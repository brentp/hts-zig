const vcf = @import("vcf.zig");
const std = @import("std");
const stdout = std.io.getStdOut().writer();

const VCF = vcf.VCF;
const InfoOrFmt = vcf.InfoOrFmt;

const fname = "tests/exac.bcf";

test "that vcf open works" {
    var ivcf = VCF.open(fname);
    defer ivcf.?.deinit();
    try std.testing.expect(ivcf != null);
}

test "that header to string works" {
    const allocator = std.testing.allocator;

    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var h = ivcf.header.tostring(allocator);
    try std.testing.expect(h != null);
    try std.testing.expect(std.mem.indexOf(u8, h.?, "#CHROM") != null);
    //try stdout.print("header:{s}\n", .{h.?});
    allocator.free(h.?);
}

test "that variant to string works" {
    const allocator = std.testing.allocator;

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
}

test "info AC int" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var fld: []const u8 = "AC";
    const allocator = std.testing.allocator;
    var ac = try variant.get(InfoOrFmt.info, i32, fld, allocator);
    try std.testing.expect(ac.len == 1);
    try std.testing.expect(ac[0] == 3);
    allocator.free(ac);
}

test "info AC float" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var fld: []const u8 = "AF";
    const allocator = std.testing.allocator;
    var af = try variant.get(InfoOrFmt.info, f32, fld, allocator);
    try std.testing.expect(af.len == 1);
    try std.testing.expect(af[0] == 6.998e-5);
    allocator.free(af);
}

test "missing INFO field gives undefined tag" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var fld: []const u8 = "MISSING_AC";
    const allocator = std.testing.allocator;
    try std.testing.expectError(vcf.HTSError.UndefinedTag, variant.get(InfoOrFmt.info, i32, fld, allocator));
}

test "format format field extraction" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    _ = ivcf.next().?;
    var variant = ivcf.next().?;
    var fld = "AD";
    const allocator = std.testing.allocator;
    var ad = try variant.get(vcf.InfoOrFmt.format, i32, fld, allocator);
    try std.testing.expect(ad.len == 8);
    try std.testing.expect(ad[0] == 7);
    try std.testing.expect(ad[2] == 2);
    try stdout.print("\nAD:{any}\n", .{ad});
    try stdout.print("\n{any}\n", .{variant});
    allocator.free(ad);
}
