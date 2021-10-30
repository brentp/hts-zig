//! This module wraps (some of) htslib's VCF parsing
//! It is very much a work in progress, but each exposed function is tested.

const std = @import("std");
const testing = std.testing;
const stdout = std.io.getStdOut().writer();
const stderr = std.io.getStdErr().writer();

const hts = @cImport({
    @cInclude("htslib_struct_access.h");
    @cInclude("htslib/hts.h");
    @cInclude("htslib/vcf.h");
    @cInclude("htslib/tbx.h");
});

pub const Field = enum {
    info,
    format,
};

/// This stores the `bcf_hdr_t` and provides convenience mthods.
pub const Header = struct {
    c: ?*hts.bcf_hdr_t,

    inline fn sync(self: Header) !void {
        const ret = hts.bcf_hdr_sync(self.c);
        if (ret != 0) {
            return ret_to_err(ret, "sync");
        }
    }
    /// add a valid string to the Header. Must contain new-line or have null
    /// terminator.
    pub fn add_string(self: Header, str: []const u8) !void {
        const ret = hts.bcf_hdr_append(self.c, &(str[0]));
        if (ret != 0) {
            return ret_to_err(ret, "header.add_string");
        }
        return self.sync();
    }

    pub fn add(self: Header, allocator: *std.mem.Allocator, fld: Field, id: []const u8, number: []const u8, typ: []const u8, description: []const u8) !void {
        // https://ziglearn.org/chapter-2/#formatting
        const h = if (fld == Field.info)
            "INFO"
        else
            "FORMAT";
        const s = try std.fmt.allocPrint(allocator, "##{s}=<ID={s},Number={s},Type={s},Description=\"{s}\">\n", .{ h, id, number, typ, description });
        defer allocator.free(s);
        return self.add_string(s);
    }

    pub fn format(self: Header, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        _ = fmt;
        _ = options;
        _ = self;
        try writer.writeAll("Header()");
    }
    pub fn tostring(self: Header, allocator: *std.mem.Allocator) ?[]u8 {
        var str = hts.kstring_t{ .s = null, .l = 0, .m = 0 };
        if (hts.bcf_hdr_format(self.c, 0, &str) != 0) {
            return null;
        }
        defer hts.free(str.s);
        defer str.s = null;
        // TODO: <strings> need copy of this, not slice.
        var sl = std.mem.sliceTo(str.s, 0);
        var result = std.mem.dupe(allocator, u8, sl) catch |err| {
            _ = err catch return "";
        };
        return result;
    }
};

/// These are the possible return values from errors in htslib calls.
pub const HTSError = error{
    IncorrectNumberOfValues,
    NotFound,
    UnexpectedType,
    UndefinedTag,
    UnknownError,
};

fn ret_to_err(
    ret: c_int,
    attr_name: []const u8,
) HTSError {
    const retval = switch (ret) {
        -10 => HTSError.IncorrectNumberOfValues,
        -3 => HTSError.NotFound,
        -2 => HTSError.UnexpectedType,
        -1 => HTSError.UndefinedTag,
        else => {
            stderr.print("[zig-hts/vcf] unknown return ({s})\n", .{attr_name}) catch {};
            return HTSError.UnknownError;
        },
    };
    return retval;
}

inline fn allele_value(val: i32) i32 {
    if (val < 0) {
        return val;
    }
    return (val >> 1) - 1;
}

pub const Allele = struct {
    val: i32,

    pub inline fn phased(a: Allele) bool {
        return (a.val & 1) == 1;
    }
    pub inline fn value(a: Allele) i32 {
        if (a.val < 0) {
            return a.val;
        }
        return (a.val >> 1) - 1;
    }
    pub fn format(self: Allele, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        _ = fmt;
        _ = options;
        const v = self.value();
        if (v < 0) {
            if (self.val == 0) {
                try writer.writeAll("./");
            } else {
                try writer.writeAll("$");
            }
        } else {
            try writer.print("{d}", .{v});
            if (self.phased()) {
                try writer.writeAll("|");
            } else {
                try writer.writeAll("/");
            }
        }
    }
};

// A genotype is a sequence of alleles.
pub const Genotype = struct {
    alleles: []i32,
    pub fn format(self: Genotype, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        for (self.alleles) |allele| {
            try (Allele{ .val = allele }).format(fmt, options, writer);
        }
    }
};

/// These are the int32 values used by htslib internally.
pub const Genotypes = struct {
    gts: []i32,
    ploidy: i32,

    pub inline fn at(self: Genotypes, i: i32) Genotype {
        var sub = self.gts[@intCast(usize, i * self.ploidy)..@intCast(usize, (i + 1) * self.ploidy)];
        return Genotype{ .alleles = sub };
    }

    pub fn format(self: Genotypes, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        _ = fmt;
        _ = options;

        try writer.writeAll("[");
        var i: i32 = 0;
        while (i * self.ploidy < self.gts.len - 1) {
            try (self.at(i)).format(fmt, options, writer);
            i += 1;
            if (i * self.ploidy < self.gts.len) {
                try writer.writeAll(", ");
            }
        }
        try writer.writeAll("]");
    }
    /// the number of alternate alleles in the genotype.
    /// for bi-allelics can think of 0 => hom-ref, 1 => het, 2=> hom-alt, -1 =>unknown
    /// ./1 == 1
    /// 0/. == 0
    /// ./. -> -1
    ///  1/1 -> 2
    ///  1/1/1 -> 3
    ///  0/2 -> 1
    ///  1/2 -> 2
    pub fn alts(self: Genotypes, allocator: *std.mem.Allocator) ![]i8 {
        var n_samples: i32 = @divTrunc(@intCast(i32, self.gts.len), self.ploidy);
        var data = try allocator.alloc(i8, @intCast(usize, n_samples));
        var i: usize = 0;
        while (i < n_samples) {
            var j: usize = 0;
            data[i] = 0;
            while (j < self.ploidy) {
                var val = std.math.min(1, std.math.max(-1, allele_value(self.gts[@intCast(usize, i * @intCast(usize, self.ploidy) + j)])));
                data[i] += @intCast(i8, val);
                j += 1;
            }
            data[i] = std.math.max(data[i], -1);
            i += 1;
        }
        return data;
    }
};

/// This provides access to the fields in a genetic variant.
pub const Variant = struct {
    c: ?*hts.bcf1_t,
    vcf: VCF,

    // deallocate memory (from htslib) for this variant. only needed if this is
    // a copy of a variant from an iterator.
    pub fn deinit(self: Variant) void {
        if (self.c != null) {
            hts.bcf_destroy(self.c);
        }
    }

    /// the 0-based start position of the variant
    pub fn start(self: Variant) i64 {
        return @as(i64, hts.variant_pos(self.c));
    }

    /// the 1-based half-open close position of the variant
    pub fn stop(self: Variant) i64 {
        return self.start() + @as(i64, hts.variant_rlen(self.c));
    }

    /// the string chromosome of the variant.
    pub fn CHROM(self: Variant) []const u8 {
        _ = hts.bcf_unpack(self.c, 4);
        var ccr = hts.bcf_hdr_id2name(self.vcf.header.c, hts.variant_rid(self.c));
        return std.mem.sliceTo(ccr, 0);
    }

    /// the reference allele
    pub fn REF(self: Variant) []const u8 {
        return std.mem.sliceTo(hts.variant_REF(self.c), 0);
    }

    /// the first alternate allele
    pub fn ALT0(self: Variant) []const u8 {
        return self.ALT(0);
    }

    /// the ith alternate allele
    pub fn ALT(self: Variant, i: i32) []const u8 {
        return std.mem.sliceTo(hts.variant_ALT(self.c, i), 0);
    }

    /// this currently returns only the first filter
    pub fn FILTER(self: Variant) []const u8 {
        _ = hts.bcf_unpack(self.c, 4);
        if (hts.variant_nflt(self.c) == 0) {
            return "PASS";
        }
        return std.mem.sliceTo(hts.variant_flt0(self.c, self.vcf.header.c), 0);
    }

    /// the variant quality
    pub fn QUAL(self: Variant) f32 {
        return hts.variant_QUAL(self.c);
    }

    /// access float or int (T of i32 or f32) in the info or format field
    /// user is responsible for freeing the returned value.
    pub fn get(self: Variant, iof: Field, comptime T: type, field_name: []const u8, allocator: *std.mem.Allocator) ![]T {
        var c_void_ptr: ?*c_void = null;

        // cfunc is bcf_get_{info,format}_values depending on `iof`.
        var cfunc = switch (iof) {
            Field.info => blk_info: {
                _ = hts.bcf_unpack(self.c, hts.BCF_UN_INFO);
                break :blk_info hts.bcf_get_info_values;
            },
            Field.format => blk_fmt: {
                _ = hts.bcf_unpack(self.c, hts.BCF_UN_FMT);
                break :blk_fmt hts.bcf_get_format_values;
            },
        };

        var n: c_int = 0;
        var typs = switch (@typeInfo(T)) {
            .ComptimeInt, .Int => .{ hts.BCF_HT_INT, i32 },
            .ComptimeFloat, .Float => .{ hts.BCF_HT_REAL, f32 },
            else => @compileError("only ints (i32, i64) and floats accepted to get()"),
        };

        var ret = cfunc(self.vcf.header.c, self.c, &(field_name[0]), &c_void_ptr, &n, typs[0]);
        if (ret < 0) {
            return ret_to_err(ret, field_name);
        }
        // typs[1] is i32 or f32
        var casted = @ptrCast([*c]u8, @alignCast(@alignOf(typs[1]), c_void_ptr));
        var data = try allocator.alloc(typs[1], @intCast(usize, n));
        @memcpy(@ptrCast([*]u8, data), casted, @intCast(usize, n * @sizeOf(typs[1])));
        return data;
    }

    /// number of samples in the variant
    pub inline fn n_samples(self: Variant) i32 {
        return hts.variant_n_samples(self.c);
    }

    /// Get the genotypes from the GT field for all samples. This allocates
    //memory that the user is expected to free.
    pub fn genotypes(self: Variant, allocator: *std.mem.Allocator) !Genotypes {
        const gts = try self.get(Field.format, i32, "GT", allocator);
        return Genotypes{ .gts = gts, .ploidy = @floatToInt(i32, @intToFloat(f32, gts.len) / @intToFloat(f32, self.n_samples())) };
    }

    pub fn format(self: Variant, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        _ = fmt;
        _ = options;
        try writer.print("Variant({s}:{d}-{d} ({s}/{s}))", .{ self.CHROM(), self.start(), self.stop(), self.REF(), self.ALT0() });
    }

    /// return a string of the variant.
    pub fn tostring(self: Variant, allocator: *std.mem.Allocator) ?[]u8 {
        var str = hts.kstring_t{ .s = null, .l = 0, .m = 0 };
        if (hts.vcf_format(self.vcf.header.c, self.c, &str) != 0) {
            return null;
        }
        var sl = std.mem.sliceTo(str.s, 0);
        defer hts.free(str.s);
        defer str.s = null;
        var result = std.mem.dupe(allocator, u8, sl) catch |err| {
            _ = err catch return "";
        };
        return result;
    }
};

/// Represents the variant (bcf or vcf) file.
/// has several convenience methods such as query and iteration.
pub const VCF = struct {
    hts: ?*hts.htsFile,
    fname: []const u8,
    header: Header,
    variant_c: ?*hts.bcf1_t,
    idx_c: ?*hts.hts_idx_t,

    /// open a vcf for reading from the given path
    pub fn open(path: []const u8) ?VCF {
        var hf = hts.hts_open(&(path[0]), "r");
        if (hf == null) {
            return null;
        }
        var h = Header{ .c = hts.bcf_hdr_read(hf.?) };
        return VCF{ .hts = hf.?, .header = h, .fname = path, .variant_c = hts.bcf_init().?, .idx_c = null };
    }

    /// set the number of decompression threads
    pub fn set_threads(self: VCF, threads: i32) void {
        hts.hts_set_threads(self.hts, @as(c_int, threads));
    }

    /// number of samples in the VCF
    pub inline fn n_samples(self: VCF) i32 {
        return hts.header_n_samples(self.header.c);
    }

    // a zig iterator over variants in the file.
    pub fn next(self: VCF) ?Variant {
        if (hts.bcf_read(self.hts, self.header.c, self.variant_c) == -1) {
            return null;
        }
        _ = hts.bcf_unpack(self.variant_c, 3);
        return Variant{ .c = self.variant_c.?, .vcf = self };
    }

    pub fn query(self: *VCF, chrom: []const u8, start: i32, stop: i32) !RegionIterator {
        if (self.idx_c == null) {
            self.idx_c = hts.hts_idx_load(&(self.fname[0]), hts.HTS_FMT_CSI);
            if (self.idx_c == null) {
                try stderr.print("[hts-zig/vcf] index not found for {any}\n", .{self.fname});
                return HTSError.NotFound;
            }
        }
        const tid = hts.bcf_hdr_name2id(self.header.c, &chrom[0]);
        if (tid == -1) {
            try stderr.print("[hts-zig/vcf] region {s} not found not found for {s}\n", .{ chrom, self.fname });
            return HTSError.NotFound;
        }

        const iter = hts.hts_itr_query(self.idx_c, tid, start, stop, hts.bcf_readrec);
        if (iter == null) {
            try stderr.print("[hts-zig/vcf] region {s}:{any}-{any} not found not found for {s}\n", .{ chrom, start + 1, stop, self.fname });
            return HTSError.NotFound;
        }
        return RegionIterator{ .itr = iter.?, .variant = Variant{ .c = self.variant_c.?, .vcf = self.* } };
    }

    /// call this to cleanup memory used by the underlying C
    pub fn deinit(self: VCF) void {
        if (self.header.c != null) {
            hts.bcf_hdr_destroy(self.header.c.?);
        }
        if (self.variant_c != null) {
            hts.bcf_destroy(self.variant_c.?);
        }
        if (self.hts != null and !std.mem.eql(u8, self.fname, "-") and !std.mem.eql(u8, self.fname, "/dev/stdin")) {
            _ = hts.hts_close(self.hts);
        }
        if (self.idx_c != null) {
            _ = hts.hts_idx_destroy(self.idx_c);
        }
    }
};

pub const RegionIterator = struct {
    itr: *hts.hts_itr_t,
    variant: Variant,

    pub fn next(self: RegionIterator) ?Variant {
        const ret = hts.hts_itr_next(hts.fp_bgzf(self.variant.vcf.hts.?), self.itr, self.variant.c, null);
        if (ret < 0) {
            hts.hts_itr_destroy(self.itr);
            return null;
        }
        // TODO: check self.variant.c.errorcode
        _ = hts.bcf_unpack(self.variant.c, 3);
        if (hts.bcf_subset_format(self.variant.vcf.header.c, self.variant.c) != 0) {
            stderr.writeAll("[hts-zig/vcf] error with bcf_subset_format\n") catch {};
            return null;
        }
        return self.variant;
    }
    // it's not necessary to call this unless iteration is stopped early.
    pub fn deinit(self: RegionIterator) void {
        _ = hts.hts_itr_destroy(self.itr);
    }
};
