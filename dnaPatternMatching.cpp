#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iomanip>

using namespace std;

const int NO_OF_CHARS = 256;
const int REPEAT_COUNT = 10000; // Number of iterations for timing

// Function to convert string to uppercase
string toUpperCase(string str) {
    string result = str;
    transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

// Display program header
void displayHeader() {
    cout << "\n===================================================" << endl;
    cout << "       DNA PATTERN MATCHING ALGORITHM ANALYZER      " << endl;
    cout << "===================================================" << endl;
}

// Display main menu
void displayMainMenu() {
    cout << "\nMAIN MENU:" << endl;
    cout << "---------" << endl;
    cout << "1. Run a single algorithm" << endl;
    cout << "2. Compare multiple algorithms" << endl;
    cout << "3. Exit program" << endl;
    cout << "\nEnter your choice (1-3): ";
}

// Print result of an algorithm
void printResult(const string& sequence, const string& pattern, const string& name, 
                const string& complexity, const vector<int>& starts, const vector<int>& ends, 
                int theoretical_time_complexity) {
    cout << "\n-------------------------------------------------" << endl;
    cout << "Algorithm: " << name << endl;
    cout << "-------------------------------------------------" << endl;
    cout << "Sequence: " << sequence << endl;
    cout << "Pattern:  " << pattern << endl;
    cout << "Number of occurrences: " << starts.size() << endl;
    
    if (!starts.empty()) {
        cout << "\nPattern occurrences:" << endl;
        cout << "-------------------" << endl;
        for (size_t i = 0; i < starts.size(); i++) {
            cout << "Match #" << (i + 1) << ":" << endl;
            cout << "  Starting index: " << starts[i] << endl;
            cout << "  Ending index:   " << ends[i] << endl;
        }
    } else {
        cout << "\nNo matches found." << endl;
    }
    
    cout << "\nAnalysis:" << endl;
    cout << "--------" << endl;
    cout << "Sequence length (n): " << setw(5) << sequence.size() << endl;
    cout << "Pattern length (m):  " << setw(5) << pattern.size() << endl;
    cout << "Time complexity:     " << complexity << endl;
    cout << "Theoretical value:   O(" << theoretical_time_complexity << ")" << endl;
    cout << "-------------------------------------------------" << endl;
}

// Validate DNA sequence
bool isValidDNA(const string& input) {
    if (input.empty()) return false;
    for (char ch : input) {
        ch = toupper(ch);
        if (ch != 'A' && ch != 'T' && ch != 'G' && ch != 'C') return false;
    }
    return true;
}

// Naive algorithm
double naive(const string& sequence, const string& pattern, bool print = true) {
    vector<int> starts, ends;
    int n = sequence.size();
    int m = pattern.size();
    
    auto start = chrono::high_resolution_clock::now();
    
    // Run algorithm multiple times to ensure measurable duration
    for (int iter = 0; iter < REPEAT_COUNT; iter++) {
        starts.clear();
        ends.clear();
        if (n == 0 || m == 0 || m > n) continue;
        
        for (int i = 0; i <= n - m; i++) {
            int j;
            for (j = 0; j < m; j++) {
                if (sequence[i + j] != pattern[j]) break;
            }
            if (j == m) {
                starts.push_back(i);
                ends.push_back(i + m - 1);
            }
        }
    }
    
    auto end = chrono::high_resolution_clock::now();
    if (print) {
        printResult(sequence, pattern, "Naive", "O(n*m)", starts, ends, n * m);
    }
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    return duration / static_cast<double>(REPEAT_COUNT); // Average time per iteration
}

// Compute LPS array for KMP
vector<int> computeLPS(const string& pattern) {
    int m = pattern.size();
    vector<int> lps(m, 0);
    int len = 0, i = 1;
    
    while (i < m) {
        if (pattern[i] == pattern[len]) {
            lps[i++] = ++len;
        } else {
            if (len != 0) len = lps[len - 1];
            else lps[i++] = 0;
        }
    }
    return lps;
}

// KMP algorithm
double kmp(const string& sequence, const string& pattern, bool print = true) {
    vector<int> starts, ends;
    int n = sequence.size();
    int m = pattern.size();
    
    auto start = chrono::high_resolution_clock::now();
    
    for (int iter = 0; iter < REPEAT_COUNT; iter++) {
        starts.clear();
        ends.clear();
        if (n == 0 || m == 0 || m > n) continue;
        
        vector<int> lps = computeLPS(pattern);
        int i = 0, j = 0;
        while (i < n) {
            if (pattern[j] == sequence[i]) {
                i++;
                j++;
            }
            if (j == m) {
                starts.push_back(i - j);
                ends.push_back(i - 1);
                j = lps[j - 1];
            } else if (i < n && pattern[j] != sequence[i]) {
                if (j != 0) j = lps[j - 1];
                else i++;
            }
        }
    }
    
    auto end = chrono::high_resolution_clock::now();
    if (print) {
        printResult(sequence, pattern, "Knuth-Morris-Pratt (KMP)", "O(n+m)", starts, ends, n + m);
    }
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    return duration / static_cast<double>(REPEAT_COUNT);
}

// Bad character heuristic for Boyer-Moore
void badCharHeuristic(const string& pattern, int size, int badchar[NO_OF_CHARS]) {
    fill_n(badchar, NO_OF_CHARS, -1);
    for (int i = 0; i < size; i++) {
        badchar[(unsigned char)pattern[i]] = i;
    }
}

// Boyer-Moore algorithm
double boyer(const string& sequence, const string& pattern, bool print = true) {
    vector<int> starts, ends;
    int n = sequence.size();
    int m = pattern.size();
    
    auto start = chrono::high_resolution_clock::now();
    
    for (int iter = 0; iter < REPEAT_COUNT; iter++) {
        starts.clear();
        ends.clear();
        if (n == 0 || m == 0 || m > n) continue;
        
        int badchar[NO_OF_CHARS];
        badCharHeuristic(pattern, m, badchar);
        int s = 0;
        
        while (s <= (n - m)) {
            int j = m - 1;
            while (j >= 0 && pattern[j] == sequence[s + j]) j--;
            if (j < 0) {
                starts.push_back(s);
                ends.push_back(s + m - 1);
                s += (s + m < n) ? m - badchar[(unsigned char)sequence[s + m]] : 1;
            } else {
                s += max(1, j - badchar[(unsigned char)sequence[s + j]]);
            }
        }
    }
    
    auto end = chrono::high_resolution_clock::now();
    if (print) {
        printResult(sequence, pattern, "Boyer-Moore", "Best: O(n/m), Worst: O(n*m)", starts, ends, n * m);
    }
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    return duration / static_cast<double>(REPEAT_COUNT);
}

// Display algorithm execution time
void displayExecutionTime(const string& algo_name, double time_ns) {
    cout << "  " << left << setw(20) << algo_name << ": ";
    
    cout << fixed << setprecision(2) << time_ns << " ms" << endl;
}

void applyOne(const string& sequence, const string& pattern) {
    cout << "\n==================================================" << endl;
    cout << "SINGLE ALGORITHM EXECUTION" << endl;
    cout << "==================================================" << endl;
    cout << "\nSelect an algorithm to run:" << endl;
    cout << "  1. Naive string matching" << endl;
    cout << "  2. Knuth-Morris-Pratt (KMP)" << endl;
    cout << "  3. Boyer-Moore" << endl;
    cout << "  4. Return to main menu" << endl;
    cout << "\nEnter your choice (1-4): ";
    
    int choice;
    cin >> choice;
    
    double execution_time = 0;
    
    switch (choice) {
        case 1:
            execution_time = naive(sequence, pattern);
            displayExecutionTime("Naive", execution_time);
            break;
        case 2:
            execution_time = kmp(sequence, pattern);
            displayExecutionTime("KMP", execution_time);
            break;
        case 3:
            execution_time = boyer(sequence, pattern);
            displayExecutionTime("Boyer-Moore", execution_time);
            break;
        case 4:
            cout << "\nReturning to main menu..." << endl;
            return;
        default:
            cout << "\nInvalid choice. Please try again." << endl;
    }
}

void applyMultiple(const string& sequence, const string& pattern) {
    cout << "\n==================================================" << endl;
    cout << "ALGORITHM COMPARISON" << endl;
    cout << "==================================================" << endl;
    cout << "\nSelect comparison option:" << endl;
    cout << "  1. Compare two algorithms" << endl;
    cout << "  2. Compare all algorithms" << endl;
    cout << "  3. Return to main menu" << endl;
    cout << "\nEnter your choice (1-3): ";

    int choice;
    cin >> choice;

    if (choice == 1) {
        cout << "\nAlgorithm options:" << endl;
        cout << "  1. Naive string matching" << endl;
        cout << "  2. Knuth-Morris-Pratt (KMP)" << endl;
        cout << "  3. Boyer-Moore" << endl;
        
        int algo1, algo2;
        cout << "\nSelect first algorithm (1-3): ";
        cin >> algo1;
        cout << "Select second algorithm (1-3): ";
        cin >> algo2;
        
        if (algo1 < 1 || algo1 > 3 || algo2 < 1 || algo2 > 3) {
            cout << "\nInvalid algorithm selection. Returning to main menu." << endl;
            return;
        }
        
        double time1 = 0, time2 = 0;
        string name1, name2;
        
        cout << "\nRunning algorithm comparison..." << endl;
        
        if (algo1 == 1) {
            time1 = naive(sequence, pattern);
            name1 = "Naive";
        } else if (algo1 == 2) {
            time1 = kmp(sequence, pattern);
            name1 = "KMP";
        } else {
            time1 = boyer(sequence, pattern);
            name1 = "Boyer-Moore";
        }
        
        if (algo2 == 1) {
            time2 = naive(sequence, pattern, false);
            name2 = "Naive";
        } else if (algo2 == 2) {
            time2 = kmp(sequence, pattern, false);
            name2 = "KMP";
        } else {
            time2 = boyer(sequence, pattern, false);
            name2 = "Boyer-Moore";
        }
        
        cout << "\n-------------------------------------------------" << endl;
        cout << "Performance Comparison Results:" << endl;
        cout << "-------------------------------------------------" << endl;
        displayExecutionTime(name1, time1);
        displayExecutionTime(name2, time2);
        cout << "-------------------------------------------------" << endl;
        
        if (time1 < time2) {
            cout << name1 << " was " << fixed << setprecision(2) 
                 << (time2 / time1) << " times faster than " << name2 << endl;
        } else if (time2 < time1) {
            cout << name2 << " was " << fixed << setprecision(2) 
                 << (time1 / time2) << " times faster than " << name1 << endl;
        } else {
            cout << "Both algorithms performed almost identically." << endl;
        }
    }else if (choice == 2) {
        cout << "\nComparing all algorithms..." << endl;

        double time_naive = naive(sequence, pattern);
        double time_kmp = kmp(sequence, pattern);
        double time_boyer = boyer(sequence, pattern);
        
        cout << "\n-------------------------------------------------" << endl;
        cout << "Performance Comparison Results:" << endl;
        cout << "-------------------------------------------------" << endl;
        displayExecutionTime("Naive", time_naive);
        displayExecutionTime("KMP", time_kmp);
        displayExecutionTime("Boyer-Moore", time_boyer);
        cout << "-------------------------------------------------" << endl;
        
        string fastest = "Naive";
        double fastest_time = time_naive;
        
        if (time_kmp < fastest_time) {
            fastest = "KMP";
            fastest_time = time_kmp;
        }
        
        if (time_boyer < fastest_time) {
            fastest = "Boyer-Moore";
            fastest_time = time_boyer;
        }
        
        cout << "The fastest algorithm was: " << fastest << endl;
    } else if (choice == 3) {
        cout << "\nReturning to main menu..." << endl;
        return;
    } else {
        cout << "\nInvalid choice. Returning to main menu." << endl;
    }
}

int main() {
    displayHeader();
    
    string sequence;
    cout << "\nEnter a DNA sequence (containing only A, T, G, C): ";
    cin >> sequence;
    
    sequence = toUpperCase(sequence);
    while (!isValidDNA(sequence)) {
        cout << "\nInvalid DNA sequence! Please use only A, T, G, C letters." << endl;
        cout << "Enter a valid DNA sequence: ";
        cin >> sequence;
        sequence = toUpperCase(sequence);
    }
    
    string pattern;
    cout << "\nEnter a DNA pattern to search for: ";
    cin >> pattern;
    
    pattern = toUpperCase(pattern);
    while (!isValidDNA(pattern)) {
        cout << "\nInvalid DNA pattern! Please use only A, T, G, C letters." << endl;
        cout << "Enter a valid DNA pattern: ";
        cin >> pattern;
        pattern = toUpperCase(pattern);
    }
    
    cout << "\nSequence: " << sequence << " (length: " << sequence.size() << ")" << endl;
    cout << "Pattern:  " << pattern << " (length: " << pattern.size() << ")" << endl;

    while (true) {
        displayMainMenu();
        
        int choice;
        cin >> choice;
        
        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "\nInvalid input. Please enter a number." << endl;
            continue;
        }
           
        switch (choice) {
            case 1:
                applyOne(sequence, pattern);
                break;
            case 2:
                applyMultiple(sequence, pattern);
                break;
            case 3:
                cout << "\nThank you for using DNA Pattern Matching Analyzer!" << endl;
                cout << "Exiting program...\n" << endl;
                return 0;
            default:
                cout << "\nInvalid choice. Please enter a number between 1 and 3." << endl;
        }
    }
    
    return 0;
}