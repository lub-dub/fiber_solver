let
  pkgs = import <nixpkgs> {};
in pkgs.mkShell {
  packages = [
    (pkgs.python3.withPackages (python-pkgs: [
      python-pkgs.ortools
      python-pkgs.black
      python-pkgs.pandas
    ]))
  ];
}
