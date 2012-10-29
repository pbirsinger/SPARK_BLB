
#!/bin/bash

mkdir sejits
cd sejits
git clone git://github.com/richardxia/asp.git

cd ~/spark
sbt/sbt assembly

cd ~
mkdir asp_extra
cd asp_extra
wget http://pypi.python.org/packages/source/c/codepy/codepy-2012.1.2.tar.gz#md5=992482e56aa3f5351a08e4c39572bc3a
tar -zxvf codepy-2012.1.2.tar.gz
cd codepy-2012.1.2
python setup.py build
python setup.py install
yum -y install numpy

cd ~
mkdir avro
cd avro
wget http://apache.osuosl.org/avro/avro-1.6.3/py/avro-1.6.3.tar.gz
tar -zxvf avro-1.6.3.tar.gz
cd avro-1.6.3
python setup.py build
python setup.py install

cd ..
mkdir java_avro
cd java_avro
wget http://www.trieuvan.com/apache/avro/avro-1.6.3/java/avro-1.6.3.jar
unzip avro-1.6.3.jar
mv org ~/avro

cd ..

wget http://jackson.codehaus.org/1.9.6/jackson-all-1.9.6.jar
unzip jackson-all-1.9.6.jar

echo "export CLASSPATH=\$CLASSPATH:.:~/avro:/root/spark/core/target/spark-core-assembly-0.4-SNAPSHOT.jar" >> ~/.bash_profile
echo "export MASTER=master@$(curl -s http://169.254.169.254/latest/meta-data/public-hostname):5050" >> ~/.bash_profile
source ~/.bash_profile

mkdir ~/models
cd ~/models
wget https://s3.amazonaws.com/halfmilEmail/comp113kmodel.avro
wget https://s3.amazonaws.com/halfmilEmail/comp250kmodel.avro
~/mesos-ec2/copy-dir /root/models

cd ~/sejits/asp/avroInter
javac -d ../avroInter/ JAvroInter.java

cp -r ~/sejits/asp/avroInter/* /root/avro
~/mesos-ec2/copy-dir /root/avro

cd ~/sejits/asp/specializers/blb
chmod +x test.sh

scalac custom_data.scala
chmod +x ~/sejits/asp/asp/jit/make_jar



