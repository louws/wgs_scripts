rm -rf /var/tmp/* /tmp/*
rm -rf /var/lib/cloud/instances/*
rm -f /var/lib/cloud/instance
rm -rf /etc/ssh/ssh_host_*
rm -f /etc/udev/rules.d/70-persistent-net.rules
find /var/log -type f -exec /bin/rm -v {} \;
touch /var/log/lastlog
cd /root
rm -f .viminfo
history -c
history -w
cat /dev/null > ~/.bash_history
sync;sync;sync
history -c
history -w
su - ec2-user -c 'rm -f .viminfo'
su - ec2-user -c 'cat /dev/null > ~/.bash_history'
sync;sync;sync
