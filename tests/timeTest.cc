#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

using namespace std;
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime()
	{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%d%m%y-%H%M%S", &tstruct);

    return buf;
	}

int main() {
	string timeNumber = currentDateTime();
    cout << "currentDateTime()=" << timeNumber << endl;
    getchar();  // wait for keyboard input
}
