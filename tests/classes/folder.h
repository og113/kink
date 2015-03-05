 /*-------------------------------------------------------------------------------------------------------------------------
 	declarations for the Folder class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __FOLDER_H_INCLUDED__
#define __FOLDER_H_INCLUDED__

#include <string>
#include <vector>
#include <utility> //for pair
#include <iostream>
#include "error.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. FilenameAttributes
	2. Errors
	3. Filename
	4. FilenameComparator
	5. Folder
	6. functions (reduceTo)
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. declarations for the FilenameAttributes class, a base class to be inherited by e.g. Filename.
		- pair<string,string> typedef
		- << pair
		- FilenameAttributes
		- <<
		- ==
-------------------------------------------------------------------------------------------------------------------------*/

typedef pair<string,string> StringPair;

template <class T, class U>
ostream& operator<<(ostream&, const pair<T,U>&);

class FilenameAttributes {
public:
	FilenameAttributes(): Directory(), Timenumber(), ID(), Suffix(), Extras() {}
	FilenameAttributes(const FilenameAttributes&);
	~FilenameAttributes() {}
	FilenameAttributes& operator=(const FilenameAttributes&);
	string 				Directory;
	string 				Timenumber;
	string 				ID;
	string 				Suffix;
	vector<StringPair> 	Extras; 
private:
	void 				copy(const FilenameAttributes&);
};

ostream& operator<<(ostream&, const FilenameAttributes&);

bool operator==(const FilenameAttributes&, const FilenameAttributes&);


/*-------------------------------------------------------------------------------------------------------------------------
	2. declarations for error structs relevant to folder
		- FilenameError
		- FilenameComparatorError
		- FolderError
-------------------------------------------------------------------------------------------------------------------------*/

class FilenameError {
public:
	class Extras: public SimpleError{
	public:
		Extras(const string& s) : Filename(s) {}		// constructor
		virtual string		message() const;			// message to be passed for printing
	private:
		string	Filename;								// filename with errors
	};
};

class FilenameComparatorError {
public:
	class LU: public SimpleError{
	public:
		LU(const string& p, const string& l, const string& u)
			 : Property(p), Lower(l), Upper(u) {}		// constructor
		virtual string		message() const;			// message to be passed for printing
	private:
		string 			Property;					// property name
		string			Lower;						// lower property
		string  		Upper;						// upper property
		
	};
};

class FolderError {
public:
	class System: public SimpleError{
	public:
		System(){}										// constructor
		virtual string		message() const;			// message to be passed for printing
	};
};

/*-------------------------------------------------------------------------------------------------------------------------
	3. declarations for the Filename class.
		- Filename
		- <<
		- >>
		- <
		- >
-------------------------------------------------------------------------------------------------------------------------*/

class Filename: public FilenameAttributes{
public:
	Filename(): FilenameAttributes() {}
	Filename(const Filename& f): FilenameAttributes(f) {}
	Filename(const string&);
	~Filename() {}
	Filename& 			operator=(const Filename&);
	Filename& 			operator=(const string&);
	string 				operator()() const;
private:
	void 				set(const string&);
};

ostream& operator<<(ostream&, const Filename&);

istream& operator>>(istream&, Filename&);

bool operator<(const Filename&, const Filename&);

bool operator>(const Filename&, const Filename&);

/*-------------------------------------------------------------------------------------------------------------------------
	4. declarations for the FilenameComparator class, which is used by the Folder class. Comparator sees the ugly details.
		- FilenameComparator
		- <<
	N.B. Upper is not used for Directory, ID or Suffix
-------------------------------------------------------------------------------------------------------------------------*/

class FilenameComparator {
public:
	FilenameComparator(): Lower(), Upper() {}
	FilenameComparator(const FilenameAttributes& l, const FilenameAttributes& u);
	FilenameComparator(const FilenameAttributes& b): Lower(b), Upper(b) {}
	FilenameComparator(const FilenameComparator&);
	~FilenameComparator() {}
	FilenameComparator& operator=(const FilenameComparator&);
	void 				set(const FilenameAttributes&, const FilenameAttributes&);
	void 				setLower(const FilenameAttributes&);
	void 				setUpper(const FilenameAttributes&);
	bool 				operator()(const Filename&) const;
	friend ostream& operator<<(ostream&, const FilenameComparator&);
private:
	FilenameAttributes 	Lower;
	FilenameAttributes 	Upper;
	bool 				check(const FilenameAttributes&,const FilenameAttributes&) const;
	void 				copy(const FilenameComparator&);
};

ostream& operator<<(ostream&, const FilenameComparator&);

/*-------------------------------------------------------------------------------------------------------------------------
	5. declarations for the Folder class.
		- Folder
		- <<
-------------------------------------------------------------------------------------------------------------------------*/

typedef vector<Filename>::iterator FolderIterator;
typedef vector<Filename>::const_iterator ConstFolderIterator;

class Folder {
public:
	Folder(const FilenameComparator&);
	Folder(const Folder&);
	Folder& operator=(const Folder&);
	~Folder() {}
	void 				set(const FilenameComparator&);
	unsigned int 		size() const;
	FolderIterator 		begin();
	ConstFolderIterator begin() const;
	FolderIterator 		end();
	ConstFolderIterator end() const;
	void 				erase(FolderIterator);
	Filename 			operator[](const int&) const;
	bool 				isPresent(const Filename&);
	void 				update();
private:
	FilenameComparator 	Comparator;
	vector<Filename> 	Filenames;
	void 				copy(const Folder&);
	void 				refresh();
	void 				sort();
};

ostream& operator<<(ostream&, const Folder&);

/*-------------------------------------------------------------------------------------------------------------------------
	6. functions acting on Filenames and Folders
		- reduceTo
-------------------------------------------------------------------------------------------------------------------------*/

void reduceTo(Folder&,const Folder&);

#endif // __FOLDER_H_INCLUDED__
