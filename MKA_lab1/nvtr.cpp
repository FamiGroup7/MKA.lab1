class nvtr {
public:
	int uzel[9], numberField;
	nvtr(int uz1, int uz2, int uz3, int uz4, int numberField_new)
	{
		uzel[0] = uz1; uzel[1] = uz2; uzel[2] = uz3; uzel[3] = uz4;
		numberField = numberField_new;
	}
	nvtr(int uz1, int uz2, int uz3, int uz4)
	{
		uzel[0] = uz1; uzel[1] = uz2; uzel[2] = uz3; uzel[3] = uz4;
		numberField = 0;
	}
	nvtr()
	{
		uzel[0] = uzel[1] = uzel[2] = uzel[3] = numberField = 0;
	}
};