#!/bin/bash

GB=$1

case $GB in
  ''|*[!0-9]*)
    echo "USAGE:"
    echo "Limit mantis memory usage to 32GB:"
    echo " $0 32 /mantis/command line"
    exit
esac

shift
MEM=$[GB*1024*1024*1024]

sudo apt install -y cgroup-tools libssl-dev

cat <<EOF > cgconfig.conf
mount {
  memory = /var/cgroups;
}
group mantis {
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

sudo umount /var/cgroups
sudo cgconfigparser -l cgconfig.conf
echo $MEM | sudo tee /var/cgroups/mantis/memory.limit_in_bytes
echo 3 | sudo tee /proc/sys/vm/drop_caches
/usr/bin/time cgexec -g memory:mantis $*
