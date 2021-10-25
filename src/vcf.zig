//! This module wraps (some of) htslib's VCF parsing
//! It is very much a work in progress, but each exposed function is tested.

const std = @import("std");
const testing = std.testing;

const hts = @cImport({
    @cInclude("htslib_struct_access.h");
    @cInclude("htslib/hts.h");
    @cInclude("htslib/vcf.h");
    @cInclude("htslib/tbx.h");
});

/// This stores the `bcf_hdr_t` and provides convenience mthods.
pub const Header = struct {
    c: *hts.bcf_hdr_t,

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

/// This provides access to the fields in a genetic variant.
pub const Variant = struct {
    c: *hts.bcf1_t,
    vcf: VCF,

    /// the 0-based start position of the variant
    pub fn start(self: Variant) i64 {
        return @as(i64, hts.variant_pos(self.c));
    }

    /// the 1-based half-open close position of the variant
    pub fn stop(self: Variant) i64 {
        return @as(i64, self.start() + @as(i64, hts.variant_rlen(self.c)));
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
        if (hts.variant_nflt(self.c) == 0) {
            return "PASS";
        }
        return std.mem.sliceTo(hts.variant_flt0(self.c, self.vcf.header.c), 0);
    }

    /// the variant quality
    pub fn QUAL(self: Variant) f32 {
        return hts.variant_QUAL(self.c);
    }

    /// access float or int (T of i32 or f32) in the info field
    /// user is responsible for freeing the returned value.
    pub fn info(self: Variant, comptime T: type, field_name: []const u8, allocator: *std.mem.Allocator) ![]T {
        _ = hts.bcf_unpack(self.c, hts.BCF_UN_INFO);
        var c_void_ptr: ?*c_void = null;

        var n: c_int = 0;
        var typs = switch (@typeInfo(T)) {
            .ComptimeInt, .Int => .{ hts.BCF_HT_INT, i32 },
            .ComptimeFloat, .Float => .{ hts.BCF_HT_REAL, f32 },
            else => @compileError("only ints (i32, i64) and floats accepted to info()"),
        };

        var ret = hts.bcf_get_info_values(self.vcf.header.c, self.c, &(field_name[0]), &c_void_ptr, &n, typs[0]);
        if (ret < 0) {
            const retval = switch (ret) {
                -10 => HTSError.IncorrectNumberOfValues,
                -3 => HTSError.NotFound,
                -2 => HTSError.UnexpectedType,
                -1 => HTSError.UndefinedTag,
                else => {
                    const stderr = std.io.getStdErr().writer();
                    try stderr.print("[zig-hts/vcf] unknown return in info({s})\n", .{field_name});
                    return HTSError.UnknownError;
                },
            };
            return retval;
        }
        // typs[1] is i32 or f32
        var casted = @ptrCast([*c]u8, @alignCast(@alignOf(typs[1]), c_void_ptr));
        var data = try allocator.alloc(typs[1], @intCast(usize, n));
        @memcpy(@ptrCast([*]u8, data), casted, @intCast(usize, n * @sizeOf(typs[1])));
        return data;
    }

    pub fn format(self: Variant, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        _ = fmt;
        _ = options;
        try writer.print("Variant({s}:{d}-{d} ({s})", .{ self.CHROM(), self.start(), self.stop(), self.REF() });
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
    hts: *hts.htsFile,
    fname: []const u8,
    header: Header,
    variant_c: *hts.bcf1_t,

    /// open a vcf for reading from the given path
    pub fn open(path: []const u8) ?VCF {
        var hf = hts.hts_open(&(path[0]), "r");
        if (hf == null) {
            return null;
        }
        var h = Header{ .c = hts.bcf_hdr_read(hf.?) };
        var v = VCF{ .hts = hf.?, .header = h, .fname = path, .variant_c = hts.bcf_init().? };

        return v;
    }

    /// set the number of decompression threads
    pub fn set_threads(self: VCF, threads: i32) void {
        hts.hts_set_threads(self.hts, @as(c_int, threads));
    }

    // a zig iterator over variants in the file.
    pub fn next(self: VCF) ?Variant {
        if (hts.bcf_read(self.hts, self.header.c, self.variant_c) == -1) {
            return null;
        }
        _ = hts.bcf_unpack(self.variant_c, 3);
        return Variant{ .c = self.variant_c, .vcf = self };
    }
};
