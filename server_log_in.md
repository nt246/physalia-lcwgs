## Logging on the AWS server from your computer
### Mac OS X and Linux users
If you are using a Mac or Linux machine, you will need to open a terminal window and then use the ```ssh``` command. ssh stands for secure shell and is a way of connecting to and interacting with remote servers. You will need to log in to the cluster using both ssh and a keyfile that has been generated for you.
Firstly, download the keyfile and open a terminal window. Then copy the keyfile into your home directory like so:
```
cp mark.pem ~
```
Then, you should be able to log in with ssh. You need to provide ssh with the path to your key, which you can do with the ```-i``` flag. This basically points to your identity file or keyfile. For example:
```
ssh -i "~/anna.pem" anna@54.245.175.86
```
Note that you will need to change the log in credentials shown here (i.e. the username and keyfile name) with your own. Also be aware that the cluster IP address will change everyday. We will update you on this each day. You might be prompted to accept an RSA key - if so, just type yes and you will log in to the cluster!

### Downloading and uploading files
Occassionally, we will need to transfer files between the cluster and our local machines. To do this, we can use a command utility called ```scp```, which stans for secure copy. It works in a similar way to ssh. Let’s try making a dummy file in our local home directory and then uploading it to our home directory on the cluster.
```
# make a file
touch test_file
# upload to cluster
scp -i "~/anna.pem" test_file mark@54.245.175.86:~/
```
Just to break this down a little we are simply copying a file, ```test_file``` in this case to the cluster. After the `:` symbol, we are specifying where on the cluster we are placing the file, here we use `~/` to specify the home directory.


Copying files back on to our local machine is just as straightforward. You can do that like so:
```
# download to local
scp -i "~/mark.pem" mark@54.245.175.86:~/test_file ./
```
Where here all we did was use scp with the cluster address and path first and the location on our computers (in our working directory) second - i.e. `./`

Alternatively, we can use Cyberduck, or similar software, to transfer data back and forth.
Open Cyberduck, and click on 'Open connection' on the top left of the screen. Select 'SFTP (SSH File Tranfer Protocol)' from the dropdown menu and copy the IP address for the day in the Server field and your username for in the Username field. Leave the Password field empty, go down to SSH Private Key and add the Private Key Carlo sent you. Press Connect and you're connected!

### Windows users
If you are using a Windows machine, you will need to log on using PuTTY since there is no native ssh client. PuTTY does not natively support the private key format (.pem) needed to login to our Amazon cloud instance. You first need to convert the private key that we gave to you to a key that PuTTY can read. PuTTY has a tool named PuTTYgen, which can convert keys to the required PuTTY format (.ppk). When you installed PuTTY, it will also have installed PuTTYgen.

First, start PuTTYgen (for example, from the Start menu, choose All Programs > PuTTY > PuTTYgen). Then select RSA and click on Load:

In the new window that pops up, Change “PuTTY Private Key Files” to “All Files” to allow you to find your pem file.

Then save your key and click on YES to dismiss the Warning as shown below.

Great, now your key file is ready and we can start Putty. In Putty, enter your user name and the IP address in the format <user_name>@<IP adress>. Make sure that 22 is given as Port and that SSH is selected.

Next, on the menu on the left, expand the “SSH” panel by clicking on the + and select “Auth”. Then, select your new putty format key file (.ppk) with Browse. Do NOT click on Open yet.

To make sure that you will not get logged out if you are inactive for a while, you can configure PuTTY to automatically send ‘keepalive’ signals at regular intervals to keep the session active. Click on Connection, then insert 180 to send keepalive signals every 3 minutes. Do NOT click on Open yet.

To avoid having to change the settings each time you log in, you can save the PuTTY session. Click on Session to get back to the basic options window. Once you are happy with all settings, give the session a name and click Save. Now you are ready to start the session with “Open”. The first time PuTTY will display a security alert dialog box that asks whether you trust the host you are connecting to. Click yes.

When you log in the next time, you can just click on the saved session and click load. If the IP address changed in the mean time (e.g. because we stopped the Amazon instance over night), you will need to replace the IP address by the new one. I would then recommend to Save the changed settings. Then simply click Open to start the session.

If the IP address did not change and you just want to login again, you can also right-click on the putty symbol in the taskbar (provided that you have pinned it to the taskbar) and select the session.

You can also tranfer files back and forth with Filezilla, a handy software to move files from a remote server such as the Amazon cloud or a cluster of your university.

Open Filezilla and choose Edit -> Settings.
Next, choose SFTP and Add the .pem key file as indicated below and click OK.
Finally, enter the IP address and the user name and when you hit enter, it should connect you. Next time, you can use the Quickconnect dropdown menu, as long as the IP address has not changed in the meantime.
Now you will see the file directory system (folders) on your local computer on the left and your folders on the amazon cloud on the right. You can now just drag and drop files from one side to the other.
