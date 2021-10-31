const vcf = @import("vcf.zig");
const std = @import("std");
const stdout = std.io.getStdOut().writer();

const VCF = vcf.VCF;
const Field = vcf.Field;

const fname = "tests/exac.bcf";
const allocator = std.testing.allocator;
const Ai32 = std.ArrayList(i32);
const Af32 = std.ArrayList(f32);

test "that vcf open works" {
    var ivcf = VCF.open(fname);
    defer ivcf.?.deinit();
    try std.testing.expect(ivcf != null);
}

test "that header to string works" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var h = ivcf.header.tostring(allocator);
    try std.testing.expect(h != null);
    try std.testing.expect(std.mem.indexOf(u8, h.?, "#CHROM") != null);
    //try stdout.print("header:{s}\n", .{h.?});
    allocator.free(h.?);
}

test "that variant to string works" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var s = variant.tostring(allocator);
    try std.testing.expect(s != null);
    try std.testing.expect(std.mem.indexOf(u8, s.?, "ENSG00000223972|") != null);
    //try stdout.print("header:{s}\n", .{s.?});
    allocator.free(s.?);
}

test "that variant attributes work" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    try std.testing.expect(std.mem.eql(u8, variant.CHROM(), "1"));
    try std.testing.expect(std.mem.eql(u8, variant.REF(), "G"));
    try std.testing.expect(std.mem.eql(u8, variant.ALT0(), "C"));
    try std.testing.expect(variant.start() == 13371);
    try std.testing.expect(variant.stop() == 13372);

    try std.testing.expect(variant.QUAL() == 608.91);
    try std.testing.expect(std.mem.eql(u8, variant.FILTER(), "PASS"));

    try std.testing.expect(std.mem.eql(u8, variant.ID(), "."));
}

test "info AC int" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var fld: []const u8 = "AC";
    var values = Ai32.init(allocator);
    defer values.deinit();

    try variant.get(Field.info, i32, &values, fld);
    try std.testing.expect(values.items.len == 1);
    try std.testing.expect(values.items[0] == 3);
}

test "info AC float" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var fld: []const u8 = "AF";
    var values = Af32.init(allocator);
    defer values.deinit();

    try variant.get(Field.info, f32, &values, fld);
    try std.testing.expect(values.items.len == 1);
    try std.testing.expect(values.items[0] == 6.998e-5);
}

test "missing INFO field gives undefined tag" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var values = Ai32.init(allocator);
    defer values.deinit();

    var fld: []const u8 = "MISSING_AC";
    try std.testing.expectError(vcf.HTSError.UndefinedTag, variant.get(Field.info, i32, &values, fld));
}

test "format format field extraction" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    _ = ivcf.next().?;
    var variant = ivcf.next().?;
    var fld = "AD";

    var ad = Ai32.init(allocator);
    defer ad.deinit();

    try variant.get(vcf.Field.format, i32, &ad, fld);
    try std.testing.expect(ad.items.len == 8);
    try std.testing.expect(ad.items[0] == 7);
    try std.testing.expect(ad.items[2] == 2);
    try stdout.print("\nAD:{any}\n", .{ad.items});
    try stdout.print("\n{any}\n", .{variant});
}

test "n_samples" {
    var gvcf = VCF.open("tests/test.snpeff.bcf").?;
    defer gvcf.deinit();
    _ = gvcf.next().?;
    var gvariant = gvcf.next().?;
    try std.testing.expect(gvariant.n_samples() == 4);
    try std.testing.expect(gvcf.n_samples() == 4);
}

test "get genotypes" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    _ = ivcf.next().?;
    var variant = ivcf.next().?;
    var gts_mem = Ai32.init(allocator);
    defer gts_mem.deinit();
    var gts = try variant.genotypes(&gts_mem);
    try stdout.print("\n{any}\n", .{gts_mem.items});
    try std.testing.expect(gts.gts.len == 8);
    try std.testing.expect(gts.ploidy == 2);
}

test "get genotypes" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    _ = ivcf.next().?;
    _ = ivcf.next().?;
    var variant = ivcf.next().?;
    var gts_mem = Ai32.init(allocator);
    defer gts_mem.deinit();

    var gts = try variant.genotypes(&gts_mem);

    try stdout.print("\ngts:{any}\n", .{gts});
    try stdout.print("\ngts.at(2):{any}\n", .{gts.at(2)});
    var alts = try gts.alts(allocator);
    try stdout.print("\nalts:{any}\n", .{alts});

    // TODO: use arraylist for alts
    allocator.free(alts);
}
test "get genotypes" {
    var ivcf = VCF.open("tests/alts.vcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var gts_mem = Ai32.init(allocator);
    defer gts_mem.deinit();
    var gts = try variant.genotypes(&gts_mem);
    var alts = try gts.alts(allocator);
    // 0/0	0/1	0/1	0/1	0/0	1/1	1/1	1/1	1/1	1/1	1/1	.	0/0	1/2	0/2	2/2
    var expected = [_]i8{ -1, 0, 1, 1, 1, 0, 2, 2, 2, 2, 2, 2, -1, 0, 2, 1, 2 };

    try std.testing.expect(std.mem.eql(i8, alts, expected[0..]));

    allocator.free(alts);
}

test "query" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    //var ivcf = VCF.open("tests/x.vcf.gz").?;
    defer ivcf.deinit();
    var chrom: []const u8 = "1";
    try std.testing.expectError(vcf.HTSError.NotFound, ivcf.query(chrom, 69269, 69270));

    chrom = "chr1";
    var iter = try ivcf.query(chrom, 69269, 69270);

    var i: i32 = 0;
    while (iter.next()) |v| {
        i += 1;
        try std.testing.expect(v.start() == 69269);
    }
    try std.testing.expect(i == 1);
}

test "add_to_header" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();

    var fld = "ASDFXX";
    var desc = "added for test";
    try ivcf.header.add(allocator, Field.info, fld, "1", "Integer", desc);
    var h = ivcf.header.tostring(allocator);
    defer allocator.free(h.?);

    try std.testing.expect(std.mem.indexOf(u8, h.?, "ASDFXX") != null);
    try std.testing.expect(std.mem.indexOf(u8, h.?, "AFXX") == null);

    _ = try ivcf.header.add(allocator, Field.format, "GGGGG", "2", "Float", desc);
}

test "set info" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();

    var v = ivcf.next().?;
    var xa = [3]i32{ 24, 48, 96 };
    var x: []i32 = xa[0..];
    var fld = "XFX";
    var desc = "added for test";
    try ivcf.header.add(allocator, Field.info, fld, "3", "Integer", desc);

    try v.set(Field.info, i32, x, fld);
    var s = v.tostring(allocator);
    try std.testing.expect(s != null);
    try std.testing.expect(std.mem.indexOf(u8, s.?, "XFX=24,48,96") != null);

    allocator.free(s.?);

    fld = "XFF";
    var xff = [_]f32{2.0};
    var xf: []f32 = xff[0..];
    try ivcf.header.add(allocator, Field.info, fld, "1", "Float", desc);

    try v.set(Field.info, f32, xf, fld);
    s = v.tostring(allocator);
    try std.testing.expect(s != null);
    try std.testing.expect(std.mem.indexOf(u8, s.?, "XFF=2") != null);
    allocator.free(s.?);
}

test "header fromstring" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    var h = ivcf.header.tostring(allocator);
    defer allocator.free(h.?);

    var hnew = vcf.Header{ .c = null };
    try hnew.from_string(h.?);
    defer hnew.deinit();

    var hs = hnew.tostring(allocator);
    defer allocator.free(hs.?);

    var l = h.?.len;

    try std.testing.expect(std.mem.eql(u8, h.?[100..l], hs.?[100..l]));
}

test "writing bcf" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    // need to add to header before copying so that the variant has the header
    // information available.
    var fld = "gnomad_popmax_af";
    try ivcf.header.add(allocator, Field.info, fld, "1", "Float", "AF from gnomad");

    var ovcf = VCF.open_mode("_o.bcf", "wb").?;
    ovcf.header.set(ivcf.header.c.?);
    ovcf.write_header();

    var xff = [_]f32{0.1};
    var xf = xff[0..];

    while (ivcf.next()) |v| {
        try v.set(Field.info, f32, xf, fld);
        try ovcf.write_variant(v);
    }

    ovcf.deinit();

    var rvcf = VCF.open("_o.bcf").?;
    defer rvcf.deinit();
    var if32 = Af32.init(allocator);
    defer if32.deinit();

    while (rvcf.next()) |*v| {
        try v.get(Field.info, f32, &if32, fld);
        try std.testing.expect(if32.items[0] == 0.1);
    }
}

test "set samples" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    var asamples = [_][]const u8{ "1094PC0009", "1094PC0012" };
    var samples = asamples[0..];
    try ivcf.set_samples(samples, allocator);
    try std.testing.expect(ivcf.n_samples() == 2);

    while (ivcf.next()) |v| {
        try std.testing.expect(v.n_samples() == 2);
    }
}
