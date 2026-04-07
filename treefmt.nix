{
  projectRootFile = "flake.nix";
  settings.global.excludes = [
    "build/**"
    "build.*/**"
    "libs/**"
    ".envrc"
  ];

  programs.nixfmt.enable = true;
  programs.meson.enable = true;

  programs.clang-format.enable = true;
  settings.formatter.clang-format.options = [ "-style=file" ];

  programs.shfmt.enable = true;
  programs.shellcheck.enable = true;

  programs.typstyle.enable = true;

  programs.typos.enable = true;
  programs.mdformat.enable = true;
  programs.mdformat.plugins = ps: [ ps.mdformat-myst ];
}
