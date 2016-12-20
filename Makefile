
all: clean
	make -C libcunbody1/
	make -C sample1/
	make -C sample2/
	make -C sample3/
	make -C sample4/


clean:
	make -C libcunbody1/ clean
	make -C sample1/ clean
	make -C sample2/ clean
	make -C sample3/ clean
	make -C sample4/ clean

