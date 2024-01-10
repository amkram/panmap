### Coding Guidelines

>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Types** ⏐ &nbsp; [globals](#global-variables
    )&nbsp;  ▪︎  &nbsp; [declaring](#variable-declaration
    )&nbsp;  ▪︎  &nbsp; [scope](#variable-scope
    )&nbsp;  ▪︎  &nbsp; [parameters](#function-parameters
)
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ⏐ &nbsp; [auto](#keyword-auto
    )&nbsp;  ▪︎  &nbsp; [integers](#integers
)

>**Organization** ⏐ &nbsp; [headers](#header-and-source-files
    )&nbsp;  ▪︎  &nbsp; [include guards](#include-guards
    )&nbsp;  ▪︎  &nbsp; [namespaces](#namespaces
)   

>&nbsp;&nbsp;&nbsp;**Formatting** ⏐ &nbsp; [brackets](#header-and-source-files
    )&nbsp;  ▪︎  &nbsp; [indentation](#indentation
    )&nbsp;  ▪︎  &nbsp; [line length](#line-length
    )&nbsp;  ▪︎  &nbsp; [case](#case-convention
)

>&nbsp;**Code quality** ⏐ &nbsp; [testing](#testing
    )&nbsp;  ▪︎  &nbsp; [documentation](#documentation
    )&nbsp;  ▪︎  &nbsp; [comments](#comments
)
---
#### Header and source files
Use header files (`.hpp`) for function declarations with definitions in a corresponding `.cpp` file. Only declare functions in headers that are exposed in the namespace. 

E.g., `pmi.hpp` declares three functions in the `pmi` namespace:

```cpp
namespace pmi {
    { ... }
    void build(seedIndex &si, const Tree *T, const size_t k, const size_t s);
    void write(std::ofstream &fout, const seedIndex &si, const Tree *T);
    void load(seedIndex &si, const Node *root, const std::ifstream &indexFile);
}
```
These three functions are implemented in `pmi.cpp` along with helper functions outside the namespace:
```cpp
#include "pmi.hpp"

/* helpers used in pmi.cpp (not in pmi namespace) */
void buildHelper(mutableTreeData &data, seedIndex &index, const Tree *T, const Node *node, const size_t k, const size_t s, const globalCoords_t &globalCoords) {
    /* ... */
}

/* ... */

/* pmi.hpp implementations, all in pmi namespace */
void pmi::build(seedIndex &index, const Tree *T, const size_t k, const size_t s) {

}
void pmi::write(std::ofstream &fout, const seedIndex &index, const Tree *T) {

}
void pmi::load(seedIndex &index, const Node *root, const std::ifstream &indexFile) {
}
```
#### Include guards
Use the `#pragma once` directive at the top of all header files. Prefer this over `#ifndef` include guards.
#### Namespaces
Encapsulate related functionality and types into namespaces as appropriate. 

Employ `using namespace foo;` directives judiciously (generally avoid).

---
#### Brackets
Use attached (same-line) braces for `if` statements, function definitions, etc.

Put a space between keywords `for`, `while`, etc. and their opening parentheses. No spaces after a function name and opening parentheses. E.g.,
```cpp

myFunction();  // good
myFunction (); // bad

for (const auto &v : container) { /* good */ }
for(const auto &v : container)  { /* bad */  }
```

Put braces on a new line after class declarations:
```cpp
class Foo() // brace on next line 
{
    // ...
}
```
#### Indentation
Indent with four spaces.
#### Line Length
No strict limit, but try to keep lines short (under ~100 characters).
#### Case convention
Use `camelCase` for names of function, variables, and data members.

Use lowercase for namespaces, preferring short names.

---

#### Global variables
Global variables should almost never be used.
#### Variable declaration
Initialize variables upon declaration when possible, prefering brace notation for initialization.
```cpp
std::vector<int> v; // Bad -- initialize on declaration
v.push_back(1);
v.push_back(2);

std::vector<int> v = {1, 2};  // Good -- v initialized (prefer brace notation)
```
#### Variable scope
Confine variables to the narrowest scope possible (declare them near their usage). E.g. in a for loop:
```cpp
for (int i = 0; i < 100; i++) {
  int num = 3; // num declared inside loop
  someFunction(num);
}
```
except when using objects, e.g.:
```cpp
Foo f; // outside loop, avoid calling constructor many times
for (int i = 0; i < 100; i++) {
    f.doSomething(); 
}
```
#### Function parameters
Pass by reference when possible, except for basic types ( int_*, char, bool, ... ).

Use `const` for all function parameters unless you need to modify them.

Place mutable parameters first in the argument list. E.g., 

```cpp
void someFunction(mutableType &result, const int32_t newValue) {
    result.value = newValue;
}
```
---
#### Keyword `auto`  
Use `auto` sparingly for local variables only. Do not use auto for file-scope or namespace-scope variables, or for class members.

Never initialize an auto-typed variable with a braced initializer list.  

Use `const auto &` in range-based loops, etc., to avoid copying. E.g.
```cpp
for (const auto &val : someVector) {
    // ...
}
```
#### Integers
`<cstdint>` defines types like `int16_t`, `uint32_t`, `int64_t`, etc. You should always use those in preference to `short`, `unsigned long long` and the like, when you need a guarantee on the size of an integer. Of the C integer types, only `int` should be used. When in doubt, `int32_t` is safe.

---

#### Testing  
We use **`boost::unit_test`** for testing. Place tests in the `src/test` directory with naming scheme `[name].test.cpp`. Place data files used in tests in `src/test/data`. The `data` directory is copied to `build` during compilation, so use relative paths like `data/[filename]` to load data in tests.

Example test file:
```cpp
#include <boost/test/unit_test.hpp>
#include "../foo.hpp"

/*
* macro that defines a new test case. Follow
* the naming convention unit_[functionName]
* for tests of a single function.
 */
BOOST_AUTO_TEST_CASE(unit_getNumber) {
    int expected = 42;
    int result = foo::getNumber(); // function being tested
    /* macro that verifies a boolean expression */
    BOOST_TEST(
        expected == result
    );
}

/*
* An "integration" test. Follow the naming convention 
* integrate_[descriptiveName] for tests that combine 
* multiple functions or logical units.
*/
BOOST_AUTO_TEST_CASE(integrate_fooLogic) {
    std::ifstream dataFile("data/testData.txt");
    std::string line;
    std::string expected = getline( dataFile, line ); 

    // Perform test logic
    std::string value = foo::getValue();
    foo::doSomething(value);
    std::string result = foo::doSomethingElse(value);

    // Verify result
    BOOST_TEST(
        expected == result
    );
}
// more test cases... 
```

#### Documentation  
Use [Doxygen](https://www.mitk.org/images/1/1c/BugSquashingSeminars$2013-07-17-DoxyReference.pdf) syntax (Javadoc style) to write comments documenting all **header files**, **public and namespace functions**, and **anything else** that is likely to be used often or by others.

Example documentation of a struct:
```cpp
/**  <- double asterisk indicates doxygen comment
 * @file 
 * @author myEmail@ucsc.edu, (Alex Kramer).
 * @author somebody@email.me, (Some Guy).

 * @brief This header file defines the Foo class.
*/
struct Foo {
public:   
    /**  <- double asterisk for auto-doc comment
     * 
     * @brief A short description of Foo.
     *
     * @details
     * Optional longer description of the struct
     * and its functionality.
     * 
     * @param [out] a an integer storing the result.
     * @param [in] s a character pointer (input parameter).
     * 
     * @return an integer with a meaningful value.
    */
    int testMe(int &a, const char *s);
};
```
Some examples of formatting Doxygen comments:
```cpp

/**
 * Use backticks for code, eg: 
 *  @brief myFunction returns an `int32_t`
 *
 * Text decoration:
 *  @brief write words in _italics_ or **bold**
 * 
 * Multi-line code blocks:
 *  @code{.cpp}
 *  std::cout << "c++ code";
 *  @endcode
*/
int32_t myFunction();
```





#### Comments

Aim for self-documenting code, with comments to explain the tricky bits.

Almost every function declaration (in headers) should have comments immediately preceding it that describe what the function does and how to use it.

Don't commit commented-out lines of code to version control.