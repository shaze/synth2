

To compile, you should be able to say `make`

If you don't have the appropriate sanitise library installed and you get an errror message like

`/opt/rh/devtoolset-9/root/usr/libexec/gcc/x86_64-redhat-linux/9/ld: cannot find -lasan`

you can just delete the `-fsanitize=address` from the Makefile

