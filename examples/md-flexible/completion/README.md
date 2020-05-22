# Command line Completions

This directory contains completion files and instructions for md-flexible.
Feel free to add completions for your favorite shell.

## zsh

1. In your `.zshrc` prepend the `zsh` subdirectory to `fpath`:
```zsh
fpath=(${PathToThisFolder}/zsh $fpath)
```

2. Initialize the auto complete system:
```zsh
autoload -U compinit && compinit
```

3. Reload your zsh
