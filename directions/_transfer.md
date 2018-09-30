## To copy data from the server

### Port connection
depending on the network you are connected you may use one of this to port to connect to the workshop server:

    -- TCP connections to the port 22 from workshop network
    -- TCP connections to the port 30022 from other networks

### linux / Mac user
 
    * Open a terminal
    * scp -r -P <PORT> <USERNAME>@gw.genome.med.kyoto-u.ac.jp:<DISTANT_PATH> <LOCAL_PATH>

### Windows

    * Open your file  transfert software (WinSCP or FilleZilla)
    * create a new connection using the following information
        .  Server : gw.genome.med.kyoto-u.ac.jp
        .  Port :   22 (from workshop network)  30022 (from the other network)
        .  USERNAME, passwd :  provided on the first day of the course
        .  Navigate to the <DISTANT_PATH>
        .  Copy the data into your <LOCAL_PATH>

