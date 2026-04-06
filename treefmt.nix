{
  projectRootFile = "flake.nix";
  settings.global.excludes = [
    "build/**"
    "build.*/**"
    "libs/**"
  ];

  programs.nixfmt.enable = true;
  programs.meson.enable = true;

  programs.clang-format.enable = true;
  settings.formatter.clang-format.options = [ "-style=file" ];
}
