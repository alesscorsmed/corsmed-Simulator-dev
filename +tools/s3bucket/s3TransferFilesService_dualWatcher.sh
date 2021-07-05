#!/bin/bash
​
# Need this to tell debconf not to use stdin as this is noninteractive script
export DEBIAN_FRONTEND=noninteractive
​
###############################################################################################
# START BASE INIT
###############################################################################################
sudo apt-get update
sudo apt-get install --force-yes --yes \
    ec2-instance-connect \
    inotify-tools \
    awscli

#echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAACAQCmpsVwX4/eouDhBmkz+uSS6CeHXZK61CFA6NWuv4HpFDRKZyrEUapK2UqirhR0Mhe3AZU0TdHDS4A1Ua+hgunOBNnG1zRTQnJjqt2RuVLDJnl9KAHE01LHIAS/4sAM8l9Y6KJSE8JAv70LcT3Uv1oL22VSmM3Qtwh+fMfbt2cCaRlUnExNNr9NOlGosFtihtuDnMeg0QwQ1r1WNb9rYJfm2RWLJBuPx9hzn18vgRKd/OmpDNtAFwa7wCOZu1s8raQU8xB9hwRuyqxaZ6dGUAx/JaYWj2Qlaj5iAwfAFYw8r39nA4KQhuLzKvsFj+fCr+ptXtCEURcQLOJLgJ1WnxRx9QXVFY1MNpRqoyxSsig2MR3LP9hQAIlyymUsDPznAFhcFA2Cs/QAYapN5Y4rUpm6jjvJdyGKXM/o42KeikXCPKipnIl5KWjROteUwbEiPjkl6s6JWGIjyUUlrm2nDBNMq27F9F1e89UupwP2bTxmA30rsiiWvBA/Xpj1M6WeR47ne/0vMKFTy4CpKtS3LFvVG97qqasQTQIthXwUAya3SECXHjaiygErbS5+wjK/uDPnnImthb9b5zJECsC4myNJaqJTAu6JTlkmVufrz79jjobZ2pzD4ItczFeDMwx39iu414cEgIwcNTp+gFOuExmiOudCS1Dp/D6DU9dPm8iR/w== corsmed-integration-key" >> /home/ubuntu/.ssh/authorized_keys

sudo cat >/usr/local/s3-upload.sh <<\EOL
#!/bin/bash
inotifywait -m /home/ubuntu/edutoolTransferToS3 -e create -e moved_to |
    while read dir action file; do
        echo "The file '$file' appeared in directory '$dir' via '$action'"
        # do something with the file
    sudo aws s3 cp /home/ubuntu/edutoolTransferToS3/$file s3://scout-images-integration.corsmed.com/ --acl public-read
    done
EOL

sudo cat >/usr/local/s3-upload-recon.sh <<\EOL
#!/bin/bash
inotifywait -m /home/ubuntu/reconimages -e create -e moved_to |
    while read dir action file; do
        echo "The file '$file' appeared in directory '$dir' via '$action'"
        # do something with the file
    sudo aws s3 cp /home/ubuntu/reconimages/$file s3://scout-images-integration.corsmed.com/ --acl public-read
    done
EOL

sudo chmod u+x /usr/local/s3-upload.sh
sudo chmod u+x /usr/local/s3-upload-recon.sh

sudo cat >/etc/systemd/system/s3-upload.service <<\EOL
[Unit]
Description=s3-upload
After=network.target
[Service]
Type=simple
Restart=on-failure
ExecStart=/bin/bash -c /usr/local/s3-upload.sh
ExecReload=/bin/kill -HUP $MAINPID
PIDFile=/var/run/s3-upload.pid
[Install]
WantedBy=multi-user.target
EOL

sudo cat >/etc/systemd/system/s3-upload-recon.service <<\EOL
[Unit]
Description=s3-upload-recon
After=network.target
[Service]
Type=simple
Restart=on-failure
ExecStart=/bin/bash -c /usr/local/s3-upload-recon.sh
ExecReload=/bin/kill -HUP $MAINPID
PIDFile=/var/run/s3-upload-recon.pid
[Install]
WantedBy=multi-user.target
EOL

sudo systemctl daemon-reload
sudo service s3-upload restart
sudo service s3-upload-recon restart
###############################################################################################