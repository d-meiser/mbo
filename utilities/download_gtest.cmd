$client = new-object System.Net.WebClient
$client.DownloadFile("http://googletest.googlecode.com/files/gtest-1.7.0.zip","gtest-1.7.0.zip")
unzip gtest-1.7.0.zip
rm gtest-1.7.0.zip
