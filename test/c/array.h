int array(void) {
	float a[9] = {0,1,2,3,4,5,6,7,8};
	float b[2][4] = {{0,1,2,3},{4,5,6,7}};
	b[0][1] = a[2];
    return 0;
}
