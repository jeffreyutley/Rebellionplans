#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

const string face[] = { "Ace", "2", "3", "4", "5", "6", "7",
                        "8", "9", "10", "Jack", "Queen", "King" }; 
const string suit[] = { "Clubs", "Diamonds", "Hearts", "Spades" };

string random_card(bool verbose=false) {
	string card;

	card = face[ rand()%13 ];
	card += " of ";
	card += suit[ rand()%4 ];

	if (verbose)
	  cout << card << "\n";

	return card;
}

int main(int argc, char *argv[])
{
	bool verbose = false;
	int seedvalue = 0;

	for (int i=1; i<argc; i++) {
	  string option = argv[i];
	  if (option.compare(0,6,"-seed=") == 0) {
	    seedvalue = atoi(&argv[i][6]);
	  } else if (option.compare("-verbose") == 0) {
	    verbose = true;
	  } else 
	    cout << "option " << argv[i] << " ignored\n";
	}

	srand(seedvalue);

	// declare a table[][] that can keep track of the
	// cards dealt for each suit -- initialize to 0
	unsigned int table[4][13] = {0};

	while (1) {

		string card = random_card(verbose);

	  // reverse engineer card suit and face (value)
	  // update accordingly: table[suit][value]++

		unsigned int i,j;
		
		for (i = 0; i < 13; ++i) {
			
			unsigned int s = face[i].size();
			
			if (card.compare(0, s, face[i]) == 0)
				break;

		}
		
		for (j = 0; j < 4; ++j) {

			unsigned int r = suit[j].size();
			unsigned int q = face[i].size();

			if (card.compare(q + 4, q + r + 4, suit[j]) == 0)
				break;

		}

		table[j][i]++;

      // break out of loop if stopping criteria met

		if (i > 9) {
			unsigned int l;
			for (l = 10; l < 13; ++l) {
				if (l == i)
					continue;

				if (table[j][l] == 0)
					break;
			}

			if (l == 13)
				break;

		}
	
	}

	// print formatted table contents to stdout 
	
	for (unsigned int a = 0; a < 4; ++a) {
		
		cout << setw(8) << right << suit[a] << " : ";
		
		for (unsigned int b = 0; b < 13; ++b)
			cout << table[a][b] << ' ';

		if ((table[a][10] != 0) && (table[a][11] != 0) && (table[a][12] != 0))
			cout << " **";

		cout << '\n';

	}

}
