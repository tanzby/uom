# uom
Unified Odometry and Mapping. A tool sets for SLAM.

## Install

Step 1: Add Bazel distribution URI as a package source

```shell
sudo apt install apt-transport-https curl gnupg -y
curl -fsSL https://bazel.build/bazel-release.pub.gpg | gpg --dearmor >bazel-archive-keyring.gpg
sudo mv bazel-archive-keyring.gpg /usr/share/keyrings
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/bazel-archive-keyring.gpg] https://storage.googleapis.com/bazel-apt stable jdk1.8" | sudo tee /etc/apt/sources.list.d/bazel.list
```


Step 2: Install and update Bazel

```shell
sudo apt update && sudo apt install bazel
```

Step 3: Install OpenGL library

```shell
sudo apt-get install libglew-dev libsdl2-dev libsdl2-image-dev libglm-dev libfreetype6-dev
```

## Use clangd for code hints

Step 1: use `ctrl + shift + P` and search `Tasks:Run Task`. Then select `Generate compile_commands.json` and start building.

Step 2:  use `ctrl + shift + P` and search `Clangd: Restart language server` and start to index the codes.
