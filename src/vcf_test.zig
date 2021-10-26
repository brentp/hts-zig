const vcf = @import("vcf.zig");
const std = @import("std");
const stdout = std.io.getStdOut().writer();

const VCF = vcf.VCF;
const InfoOrFmt = vcf.InfoOrFmt;

const fname = "tests/exac.bcf";
const allocator = std.testing.allocator;

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
}

test "info AC int" {
    var ivcf = VCF.open("tests/exac.bcf").?;
    defer ivcf.deinit();
    var variant = ivcf.next().?;
    var fld: []const u8 = "AC";
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
    try std.testing.expectError(vcf.HTSError.UndefinedTag, variant.get(InfoOrFmt.info, i32, fld, allocator));
}

test "format format field extraction" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    _ = ivcf.next().?;
    var variant = ivcf.next().?;
    var fld = "AD";
    var ad = try variant.get(vcf.InfoOrFmt.format, i32, fld, allocator);
    try std.testing.expect(ad.len == 8);
    try std.testing.expect(ad[0] == 7);
    try std.testing.expect(ad[2] == 2);
    try stdout.print("\nAD:{any}\n", .{ad});
    try stdout.print("\n{any}\n", .{variant});
    allocator.free(ad);
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
    var gts = try variant.genotypes(allocator);
    try std.testing.expect(gts.gts.len == 8);
    try std.testing.expect(gts.ploidy == 2);
    allocator.free(gts.gts);
}

test "get genotypes" {
    var ivcf = VCF.open("tests/test.snpeff.bcf").?;
    defer ivcf.deinit();
    _ = ivcf.next().?;
    _ = ivcf.next().?;
    var variant = ivcf.next().?;
    var gts = try variant.genotypes(allocator);

    try stdout.print("\ngts:{any}\n", .{gts});
    try stdout.print("\ngts.at(2):{any}\n", .{gts.at(2)});

    allocator.free(gts.gts);
}
