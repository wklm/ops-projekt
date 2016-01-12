// Implementieren Sie den Simplex Algorithmus. Schreiben Sie dafur eine Funkti-
//on mit der Signatur
//double* lpsolve(int n, double *c, int k, double **A, double *b).
//39
//n gibt die Anzahl der Variablen an, k die Anzahl der Nebenbedingung. Die
//Vektoren c, b und die Matrix A sind entsprechend der Vorlesung zu ubergeben,
//sodass das Lineare Programm in Standardform
//c
//T
//x ? max!
//    Ax ? b
//x ? 0
//ist.
//    Ruckgabe wird ein Array mit den entsprechenden Loesungswerten sein
// in der ersten Zeile stehen zwei mit Abstand getrennte Zahlen n, k, die wie
//oben zu verstehen sind.
// in der zweiten Zeile stehen n Fliesskommazahlen, die
// mit Absanden getrennt sind, und den Vektor c reprasentieren.
//nun folgen k Zeilen mit jeweils n + 1 Zahlen. Die ersten n Zahlen bilden
//eine Zeile der Matrix A, die letzte steht fur den zugehoerigen Eintrag im
//Vektor b (Jede Zeile entspricht daher einer Nebenbedingung.)
//Beachten Sie zum Einlesen eines ahnlichen Files auch das Beispielprogramm
//im Anhang.

if (typeof console === "undefined"){
    console={};
    console.log = function(){
        return;
    }
}

var data="2 3 \n"+
    "300 500 \n"+
    "1	2	170 \n"+
    "1	1	150 \n"+
    "0	2	180\n";

var data2="2 3 \n"+
    "3 5 \n"+
    "1	0	4 \n"+
    "0	2	12 \n"+
    "3	2	18 \n";

//function double * lpsolve(int n, double * c, int k, double ** A, double * b) {
function parseSolve(data, steps)
{
    var lines = data.split('\n');

    //first line: n,k
    var first=lines[0].trim().split(/\s+/);
    var n=parseInt(first[0],10);
    var k=parseInt(first[1],10);

    //second line: n float values
    var cValsText=lines[1].trim().split(/\s+/);
    var cVals =[];
    for(i=0;i<cValsText.length;i++)
        cVals.push(-parseFloat(cValsText[i],10));

    //remaining lines
    var b=new Array();
    var A = create2DArray(k);
    for(var i = 0; i < k; i++){
        var line=lines[i+2];
        var tokens=line.trim().split(/\s+/);
        for(var j=0;j<tokens.length-1;j++)
            A[i][j]=parseFloat(tokens[j]);
        b[i]=parseFloat(tokens[tokens.length-1]);
    }

    printVector("c", cVals);
    printMatrix("A", A);
    printVector("b", b);

    for(i=0;i<k;i++)
    {
        for(j=0;j<k;j++)
        {
            if(i==j) A[i][j+n]=1;
            else A[i][j+n]=0;

        }
    }

    printVector("c", cVals);
    printMatrix("A", A);
    printVector("b", b);
    return lpsolve(n, cVals, k, A, b, steps);
}

function printVector(name, v)
{
    console.log("Vector "+name +": ");
    var line="";
    for(var i = 0; i < v.length; ++i)
        line+=v[i] + "\t";
    console.log(line);
}

function printMatrix(name, m)
{
    console.log("Matrix "+name +": ");
    for(var i = 0; i < m.length; ++i) {
        var line="";
        for(var j = 0; j < m[i].length; ++j)
            line+=m[i][j] + "\t";
        console.log(line);
    }
}

function Result()  {
    this.step=undefined;
    this.A=undefined;
    this.c=undefined;
    this.b=undefined;
    this.pivotCol=undefined;
    this.pivotRow=undefined;
    this.bv=undefined;
}


//passing arrays and variables by reference, so they can be changed if needed
function dualSimplex(n, c, k, A, b) {
    //check if all bi in array b are negative (dual simplex required if so)
    var negCounter = 0;
    for(var i = 0; i < k; ++i) {
        if(b[i] < 0) {
            ++negCounter;
        }
    }
    if(negCounter == k) {
        //swap the values of n and k (n is actually n + k)
        var tempK = k;
        k = n - k;
        n = tempK + k;
        //keep old b values in a temp array
        var tempB = b;
        //update b (corresponds to the dual c matrix from the OPS script)
        b = new double[k + 1]();
        for(var i = 0; i < k; ++i) {
            b[i] = c[i];
        }
        //update c (corresponds to the dual b matrix from the OPS script)
        c = new double[n]();
        for(var i = 0; i < tempK; ++i) {
            c[i] = tempB[i];
        }
        //transpone matrix A
        var tempA = A;
        var A = create2DArray(k)
        for(var i = 0; i < k; ++i) {
            var j = 0;
            while(j < n - k) {
                A[i][j] = -tempA[j][i];
                ++j;
            }
            A[i][j + i] = 1;
        }
        console.log("**DUAL SIMPLEX**");
        //print(n, c, k, A, b);
        return true;
    }
    return false;
}

function  lpsolve(n, c, k, A, b, steps) {
    var dual = dualSimplex(n, c, k, A, b);
    console.log(dual);

    var min;
    var solution=new Array(); //n+1; //double * solution = new double [n + 1];
    var bv=new Array(); //k //int * rowVariableTracker = new int [k];
    for(i=n;i<n+k;i++)
        c[i]=0;
    b[k]=0;

    for(var i = 0; i < k; i++)
        bv[i] = n + i;

    var curStep=0;
    while(curStep<steps&&(min = findMin(c, n+k)) < 0) {
        curStep++;

        var pivotColumn;
        for(var i = 0; i < n+k; ++i)
            if(c[i] == min)
                pivotColumn = i;

        var pivotRowSet = false;
        var pivotRowDiv;
        var pivotRow;
        var pivot;
        for(var i = 0; i < k; ++i) {
            if(A[i][pivotColumn] != 0 && A[i][pivotColumn] > 0) {
                if(!pivotRowSet) {
                    pivotRowDiv = b[i] / A[i][pivotColumn];
                    pivotRow = i;
                    pivot = A[i][pivotColumn];
                    pivotRowSet = true;
                }
                if(b[i] / A[i][pivotColumn] < pivotRowDiv) {
                    pivotRowDiv = b[i] / A[i][pivotColumn];
                    pivotRow = i;
                    pivot = A[i][pivotColumn];
                }
            }
        }
        console.log("-> Pivot Column:" + pivotColumn);
        console.log("-> Pivot Row:   " + pivotRow );

        bv[pivotRow] = pivotColumn;
        var divisor;
        var functionDivisor = c[pivotColumn] / pivot;

        for(var j = 0; j < k; ++j) {
            if(j == pivotRow) {
                continue;
            } else {
                divisor = A[j][pivotColumn] / pivot;
            }
            for(var f = 0; f < n+k; ++f) {
                A[j][f] = A[j][f] - (divisor * A[pivotRow][f]);
            }
            b[j] = b[j] - (divisor * b[pivotRow]);
        }

        b[k]=b[k]-(functionDivisor * b[pivotRow]);
        //b[0] = b[0] - (functionDivisor * b[pivotRow]);
        for(var p = 0; p < n+k; ++p) {
            c[p] = c[p] - (functionDivisor * A[pivotRow][p]);
        }
        b[pivotRow] = b[pivotRow] / pivot;
        for(var s = 0; s < n+k; ++s) {
            A[pivotRow][s] = A[pivotRow][s] / pivot;
        }
        printVector("C",c);
        printMatrix("A",A);
        printVector("B",b);
    }
    r=new Result();
    r.A=A;
    r.c=c;
    r.b=b;
    r.pivotCol=pivotColumn;
    r.pivotRow=pivotRow;
    r.step=step;
    r.bv=bv;
    printVector("C",c);
    printMatrix("A",A);
    printVector("B",b);
    return r;
}

//parse(data2);

//http://stackoverflow.com/questions/966225/how-can-i-create-a-two-dimensional-array-in-javascript
function create2DArray(rows) {
    var arr = [];

    for (var i=0;i<rows;i++) {
        arr[i] = [];
    }

    return arr;
}

function findMin(c, n) {
    var temp = c[0];
    for(var i = 0; i < n; ++i)
        if(c[i] < temp)
            temp = c[i];
    return temp;
}