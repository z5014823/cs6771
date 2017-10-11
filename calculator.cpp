#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stack>
#include <utility>
#include <queue>
#include <vector>

// TODO, how to destruct stacks, queues, vectors containing objects.

void operate(std::stack<std::pair<double, std::string> >&, std::string&);

int main(int argc, const char * argv[]) {
    // Set double to 3 decimal places.
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(3);
    
    // Open input .txt file.
    std::ifstream in;
    in.open(argv[1]);
    
    // Read input .txt file.
    std::string s;
    std::stack<std::pair<double, std::string> > calc;
    while (in >> s) {
        // If token is a number, push to stack.
        if (isdigit(s[0])) {
            // If number is a double.
            if (s.find('.') != std::string::npos) {
                calc.push(std::make_pair(std::stoi(s), "double"));
                std::cout << "Pushed " << (calc.top()).first << " of type " << (calc.top()).second << "." << std::endl;
            }
            // If number is an integer.
            else {
                calc.push(std::make_pair(std::stoi(s), "int"));
                std::cout << "Pushed " << (int)(calc.top()).first << " of type " << (calc.top()).second << "." << std::endl;
            }
        }
        // If token is "repeat".
        else if (s == "repeat") { // "repeat" token is not part of main operations, otherwise circular call to outside nest needed.
            int n = (int)((calc.top()).first);
            
            // If top of stack is the number 0, can't repeat 0 times. DO NOT CONTINUE. Must process all tokens before "endrepeat".
            if (n == 0) std::cout << "Top of stack is the number 0, can't repeat 0 times." << std::endl;
            
            // Store all tokens preceding "endrepeat" in a vector of strings.
            std::vector<std::string> tokens;
            while (in >> s && s != "endrepeat") tokens.push_back(s);
            
            // If no tokens were added to the vector, nothing to repeat.
            if (tokens.empty()) {
                //~tokens();
                std::cout << "No tokens were added to the vector, nothing to repeat." << std::endl;
                continue;
            } else {
                // Perform n repetitions.
                for (int i = 0; i < n; ++i) {
                    // Iterate over stored tokens.
                    for (auto s : tokens) {
                        // If token is a number, push to stack. We can use the same stack as in other nests. Consider if 0 repetitions are made, then no elements in the tokens vector will be processed.
                        if (isdigit(s[0])) {
                            // If number is a double.
                            if (s.find('.') != std::string::npos) {
                                calc.push(std::make_pair(std::stoi(s), "double"));
                                std::cout << "Pushed " << (calc.top()).first << " of type " << (calc.top()).second << "." << std::endl;
                            }
                            // If number is an integer.
                            else {
                                calc.push(std::make_pair(std::stoi(s), "int"));
                                std::cout << "Pushed " << (int)(calc.top()).first << " of type " << (calc.top()).second << "." << std::endl;
                            }
                        }
                        // If token is "add", "sub", "mult", "div", "sqrt", "pop" or "reverse".
                        else {
                            operate(calc, s);
                        }
                    }
                }
                //~tokens();
            }
        // If token is "add", "sub", "mult", "div", "sqrt", "pop" or "reverse".
        } else {
            operate(calc, s);
        }
    }
    in.close();
    
    return EXIT_SUCCESS;
}

void operate(std::stack<std::pair<double, std::string> >& calc, std::string& s) {
    if (s == "add") {
        // Check that the first element is available.
        if (calc.empty()) {
            std::cout << s << " requires 2 elements. 0 available." << std::endl;
        } else {
            double a = (calc.top()).first;
            std::string atype = (calc.top()).second;
            calc.pop();
            
            // Check that the second element is available.
            if (calc.empty()) {
                calc.push(make_pair(a, atype));
                std::cout << s << " requires 2 elements. 1 available." << std::endl;
            } else {
                double b = (calc.top()).first;
                std::string btype = (calc.top()).second;
                calc.pop();
                
                if (atype == "int" && btype == "int") {
                    std::cout << (int)a << " + " << (int)b << " = " << (int)(a + b) << std::endl;
                } else if (atype == "int") {
                    std::cout << (int)a << " + " << b << " = " << a + b << std::endl;
                } else if (btype == "int") {
                    std::cout << a << " + " << (int)b << " = " << a + b << std::endl;
                } else {
                    std::cout << a << " + " << b << " = " << a + b << std::endl;
                }
            }
        }
    } else if (s == "sub") {
        // Check that the first element is available.
        if (calc.empty()) {
            std::cout << s << " requires 2 elements. 0 available." << std::endl;
        } else {
            double a = (calc.top()).first;
            std::string atype = (calc.top()).second;
            calc.pop();
            
            // Check that the second element is available.
            if (calc.empty()) {
                calc.push(make_pair(a, atype));
                std::cout << s << " requires 2 elements. 1 available." << std::endl;
            } else {
                double b = (calc.top()).first;
                std::string btype = (calc.top()).second;
                calc.pop();
                
                if (atype == "int" && btype == "int") {
                    std::cout << (int)a << " - " << (int)b << " = " << (int)(a - b) << std::endl;
                } else if (atype == "int") {
                    std::cout << (int)a << " - " << b << " = " << a - b << std::endl;
                } else if (btype == "int") {
                    std::cout << a << " - " << (int)b << " = " << a - b << std::endl;
                } else {
                    std::cout << a << " - " << b << " = " << a - b << std::endl;
                }
            }
        }
    } else if (s == "mult") {
        // Check that the first element is available.
        if (calc.empty()) {
            std::cout << s << " requires 2 elements. 0 available." << std::endl;
        } else {
            double a = (calc.top()).first;
            std::string atype = (calc.top()).second;
            calc.pop();
            
            // Check that the second element is available.
            if (calc.empty()) {
                calc.push(make_pair(a, atype));
                std::cout << s << " requires 2 elements. 1 available." << std::endl;
            } else {
                double b = (calc.top()).first;
                std::string btype = (calc.top()).second;
                calc.pop();
                
                if (atype == "int" && btype == "int") {
                    std::cout << (int)a << " * " << (int)b << " = " << (int)(a * b) << std::endl;
                } else if (atype == "int") {
                    std::cout << (int)a << " * " << b << " = " << a * b << std::endl;
                } else if (btype == "int") {
                    std::cout << a << " * " << (int)b << " = " << a * b << std::endl;
                } else {
                    std::cout << a << " * " << b << " = " << a * b << std::endl;
                }
            }
        }
    } else if (s == "div") {
        // Check that the first element is available.
        if (calc.empty()) {
            std::cout << s << " requires 2 elements. 0 available." << std::endl;
        } else {
            double a = (calc.top()).first;
            std::string atype = (calc.top()).second;
            calc.pop();
            
            // Check that the second element is available.
            if (calc.empty()) {
                calc.push(make_pair(a, atype));
                std::cout << s << " requires 2 elements. 1 available." << std::endl;
            } else {
                double b = (calc.top()).first;
                std::string btype = (calc.top()).second;
                calc.pop();
                
                if (atype == "int" && btype == "int") {
                    if ((int)a % (int)b == 0) {
                        std::cout << (int)a << " / " << (int)b << " = " << (int)(a / b) << std::endl;
                    } else {
                        std::cout << (int)a << " / " << (int)b << " = " << (a + 0.0) / b << std::endl;
                    }
                } else if (atype == "int") {
                    std::cout << (int)a << " - " << b << " = " << a / b << std::endl;
                } else if (btype == "int") {
                    std::cout << a << " / " << (int)b << " = " << a / b << std::endl;
                } else {
                    std::cout << a << " / " << b << " = " << a / b << std::endl;
                }
            }
        }
    } else if (s == "sqrt") {
        double a = (calc.top()).first;
        std::string atype = (calc.top()).second;
        calc.pop();
        if (atype == "int" && sqrt(a) - (int)sqrt(a) == 0) {
            std::cout << "sqrt " << (int)a << " = " << (int)sqrt(a) << std::endl;
        } else {
            std::cout << "sqrt " << a << " = " << sqrt(a) << std::endl;
        }
    } else if (s == "pop") {
        calc.pop();
    } else if (s == "reverse") {
        int n = (int)((calc.top()).first);
        std::cout << "Reversing the last " << n << " elements of the stack." << std::endl;
        std::queue<std::pair<double, std::string> > rev;
        int i = 0;
        while (i++ < n && !calc.empty()) {
            rev.push(calc.top());
            calc.pop();
        }
        while (!rev.empty()) {
            calc.push(rev.front());
            rev.pop();
        }
        //~rev();
    }
    // Additional print token for debugging.
    else if (s == "print") {
        std::stack<std::pair<double, std::string> > tmp;
        while (!calc.empty()) {
            tmp.push(calc.top());
            calc.pop();
        }
        std::cout << "Stack elements: ";
        while (!tmp.empty()) {
            std::cout << (tmp.top()).first << " of type " << (tmp.top()).second << ". ";
            calc.push(tmp.top());
            tmp.pop();
        }
        std::cout << std::endl;
        //~tmp();
    }
    // Handling bad input. Currently ineffective if first character of bad input is a number.
    else {
        std::cout << "Ignoring bad input." << std::endl;
    }
}
