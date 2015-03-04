 /*
 	declarations for the Folder class and dependencies
 */
 
#ifndef __FOLDER_H_INCLUDED__
#define __FOLDER_H_INCLUDED__

#include <string>
#include <vector>
#include <utility> //for pair
#include <iostream>

using namespace std;

/*
	declarations for the FilenameAttributes class, a base class to be inherited by e.g. Filename.
		- pair<string,string> typedef
		- FilenameAttributes
		- <<
*/

typedef pair<string,string> StringPair;

class FilenameAttributes {
public:
	FilenameAttributes(): Directory(), Timenumber(), ID(), Suffix(), Extras() {}
	FilenameAttributes(const FilenameAttributes&);
	FilenameAttributes& operator=(const FilenameAttributes&);
	~FilenameAttributes() {}
	string Directory;
	string Timenumber;
	string ID;
	string Suffix;
	vector<StringPair> Extras; 
private:
	void copy(const FilenameAttributes&);
};

ostream& operator<<(ostream&, const FilenameAttributes&);

/*
	declarations for the Filename class.
		- Filename
		- <<
*/

class Filename: public FilenameAttributes{
public:
	Filename(): FilenameAttributes() {}
	Filename(const Filename& f): FilenameAttributes(f) {}
	Filename(const string&);
	Filename& operator=(const Filename&);
	~Filename() {}
	string operator()() const;
};

ostream& operator<<(ostream&, const Filename&);

/*
	declarations for the FilenameComparator class, which is used by the Folder class. Comparator sees the ugly details.
		- FilenameComparator
		- <<
	N.B. Upper is not used for Directory, ID or Suffix
*/

class FilenameComparator {
public:
	FilenameComparator(): Lower(), Upper() {}
	FilenameComparator(const FilenameAttributes& l, const FilenameAttributes& u);
	FilenameComparator(const FilenameAttributes& b): Lower(b), Upper(b) {}
	FilenameComparator(const FilenameComparator&);
	FilenameComparator& operator=(const FilenameComparator&);
	void set(const FilenameAttributes&, const FilenameAttributes&);
	void setLower(const FilenameAttributes&);
	void setUpper(const FilenameAttributes&);
	bool operator()(const Filename&) const;
	~FilenameComparator() {}
	friend ostream& operator<<(ostream&, const FilenameComparator&);
private:
	FilenameAttributes Lower;
	FilenameAttributes Upper;
	bool check(const FilenameComparator&,const FilenameComparator&);
	void copy(const FilenameComparator&);
};

ostream& operator<<(ostream&, const FilenameComparator&);

/*
	declarations for the Folder class.
		- Folder
		- <<
*/

class Folder {
public:
	Folder(const FilenameComparator&);
	Folder(const Folder&);
	Folder& operator=(const Folder&);
	void set(const FilenameComparator&);
	int size();
	void update();
	Filename& operator[](const int&);
	bool isPresent(const Filename&);
private:
	FilenameComparator Comparator;
	vector<Filename> Filenames;
	void copy(const Folder&);
	void refresh();
}

ostream& operator<<(ostream&, const Folder&);

#endif // __FOLDER_H_INCLUDED__
