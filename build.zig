const Builder = @import("std").build.Builder;

pub fn build(b: *Builder) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.

    // Standard release options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall.
    const mode = b.standardReleaseOptions();
    // https://ziglang.org/documentation/0.8.0/#Building-a-Library

    const lib = b.addStaticLibrary("zigvcf", "src/vcf.zig");
    lib.linkLibC();
    lib.addIncludeDir("src/");
    lib.addCSourceFile("src/htslib_struct_access.c", &[_][]const u8{});
    lib.linkSystemLibrary("hts");
    lib.linkSystemLibrary("z");
    lib.setBuildMode(mode);
    lib.install();

    var tests = b.addTest("src/vcf_test.zig");
    tests.setBuildMode(mode);
    tests.linkLibC();
    tests.linkSystemLibrary("hts");
    tests.addIncludeDir("src/");
    tests.addCSourceFile("src/htslib_struct_access.c", &.{});

    const test_step = b.step("test", "Run the tests");
    test_step.dependOn(&tests.step);
}
