#!/bin/bash

# The example script for building mantis on a limited memory of 4GBs assuming we are in the root of <mantis_dir> is:
# bash script/limitMem.sh 4 build/src/mantis build -s 25 -i <input_squeakr_list_file> -o <mantis index>


GB=$1

case $GB in
  ''|*[!0-9]*)
    echo "USAGE:"
    echo "Limit a program's memory usage to <n> GB:"
    echo " $0 <n> <program command line>"
    exit
esac

shift
MEM=$[GB*1024*1024*1024]

CGCONFIGPARSER=`which cgconfigparser`
if [ ! -n "$CGCONFIGPARSER" ]; then
    for guess in /sbin /usr/sbin /usr/local/sbin; do
        if [ ! -x $guess/cgconfigparser ]; then
            CGCONFIGPARSER=$guess/cgconfigparser
            break
        fi
    done
fi
if [ ! -n "$CGCONFIGPARSER" ]; then
    echo "This script requires cgconfigparser."
    echo "Try installing the cgroup-tools package."
    echo "On Debian-based systems, such as Ubuntu:"
    echo "  sudo apt install -y cgroup-tools libssl-dev"
    exit 1
fi

CGEXEC=`which cgexec`
if [ ! -n "$CGEXEC" ]; then
    for guess in /sbin /usr/sbin /usr/local/sbin; do
        if [ ! -x $guess/cgexec ]; then
            CGEXEC=$guess/cgexec
            break
        fi
    done
fi
if [ ! -n "$CGEXEC" ]; then
    echo "This script requires cgexec."
    echo "Try installing the cgroup-tools package."
    echo "On Debian-based systems, such as Ubuntu:"
    echo "  sudo apt install -y cgroup-tools libssl-dev"
    exit 1
fi

GROUPNAME=`mktemp -u cgroupname.XXXXXX`
if [ $? -ne 0 ]; then
    echo "Failed to generate ephemeral cgroup name"
    exit 1
fi

CONFIGFILE=`mktemp -t cgconfig.XXXXXX`
if [ $? -ne 0 ]; then
    echo "Failed to create config file."
    exit 1
fi

MOUNTDIR=`mktemp -t -d cgroup.XXXXXX`
if [ $? -ne 0 ]; then
    echo "Failed to create mount directory."
    exit 1
fi

cat <<EOF > $CONFIGFILE
mount {
  memory = $MOUNTDIR;
}
group $GROUPNAME {
  perm {
    task {
      uid = `id -u`;
      gid = `id -g`;
      fperm = 770;
    }
    admin {
      uid = `id -u`;
      gid = `id -g`;
      fperm = 770;
    }
  }
  memory {
  }
}
EOF

if [ $? -ne 0 ]; then
    echo "Failed to write config file"
    exit 1
fi

#sudo umount /var/cgroups
# No error check here.

sudo $CGCONFIGPARSER -l $CONFIGFILE
if [ $? -ne 0 ]; then
    echo "Failed to setup the mantis cgroup"
    exit 1
fi

echo $MEM | sudo tee $MOUNTDIR/$GROUPNAME/memory.limit_in_bytes
MEMCHECK=`cat $MOUNTDIR/$GROUPNAME/memory.limit_in_bytes`
if [ $MEM != $MEMCHECK ]; then
    echo "Failed to set memory limit."
    sudo umount $MOUNTDIR
    exit 1
fi

echo 3 | sudo tee /proc/sys/vm/drop_caches

/usr/bin/time $CGEXEC -g memory:$GROUPNAME $*
result=$?

sudo umount $MOUNTDIR
if [ $? -ne 0 ]; then
    echo "Warning: failed to cleanup/unmount $MOUNTDIR"
fi

rmdir $MOUNTDIR
if [ $? -ne 0 ]; then
    echo "Warning: failed to delete $MOUNTDIR"
fi

exit $result
