// Implementieren Sie den Simplex Algorithmus. Schreiben Sie dafur eine Funkti- ®
//on mit der Signatur
//double* lpsolve(int n, double *c, int k, double **A, double *b).
//39
//n gibt die Anzahl der Variablen an, k die Anzahl der Nebenbedingung. Die
//Vektoren c, b und die Matrix A sind entsprechend der Vorlesung zu ubergeben, ®
//sodass das Lineare Programm in Standardform
//c
//T
//x ? max!
//    Ax ? b
//x ? 0
//ist.
//    Ruckgabe wird ein Array mit den entsprechenden L ® osungswerten sein
//ñ in der ersten Zeile stehen zwei mit Abstand getrennte Zahlen n, k, die wie
//oben zu verstehen sind.
//ñ in der zweiten Zeile stehen n Flieﬂkommazahlen, die mit Absanden ge- ®
//trennt sind, und den Vektor c reprasentieren. ®
//ñ nun folgen k Zeilen mit jeweils n + 1 Zahlen. Die ersten n Zahlen bilden
//eine Zeile der Matrix A, die letzte steht fur den zugeh ® origen Eintrag im ®
//Vektor b (Jede Zeile entspricht daher einer Nebenbedingung.)
//ñ Beachten Sie zum Einlesen eines ahnlichen Files auch das Beispielprogramm ®
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

    b.splice(0,0,0); //TODO: ?
    printVector("c", cVals);
    printMatrix("A", A);
    printVector("b", b);
    return lpsolve(n+k, cVals, k, A, b, steps);
    printVector("c", cVals);
    printMatrix("A", A);
    printVector("b", b);
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

var Result ={
    step: undefined,
    A: undefined,
    c:undefined,
    b:undefined,
    pivotCol:undefined,
    pivotRow:undefined
}

function  lpsolve(n, c, k, A, b, steps) {


    var min;
    var solution=new Array(); //n+1; //double * solution = new double [n + 1];
    var rowVariableTracker=new Array(); //k //int * rowVariableTracker = new int [k];

    for(var i = 0; i < k; ++i) {
        rowVariableTracker[i] = n - k + i + 1;
    }
    var curStep=0;
    while(curStep<steps&&(min = findMin(c, n)) < 0) {
        curStep++;
        var pivotColumn;
        for(var i = 0; i < n; ++i)
            if(c[i] == min)
                pivotColumn = i;


        var pivotRowSet = false;
        var pivotRowDiv;
        var pivotRow;
        var pivot;
        for(var i = 1; i < k + 1; ++i) {
            if(A[i - 1][pivotColumn] != 0) {
                if(!pivotRowSet) {
                    pivotRowDiv = b[i] / A[i - 1][pivotColumn];
                    pivotRow = i - 1;
                    pivot = A[i - 1][pivotColumn];
                    pivotRowSet = true;
                }
                if(b[i] / A[i - 1][pivotColumn] < pivotRowDiv) {
                    pivotRowDiv = b[i] / A[i - 1][pivotColumn];
                    pivotRow = i - 1;
                    pivot = A[i - 1][pivotColumn];
                }
            }
        }
        console.log("-> Pivot Column:" + pivotColumn);
        console.log("-> Pivot Row:   " + pivotRow );

        rowVariableTracker[pivotRow] = pivotColumn + 1;
        var divisor;
        var functionDivisor = c[pivotColumn] / pivot;

        for(var j = 0; j < k; ++j) {
            if(j == pivotRow) {
                continue;
            } else {
                divisor = A[j][pivotColumn] / pivot;
            }
            for(var f = 0; f < n; ++f) {
                A[j][f] = A[j][f] - (divisor * A[pivotRow][f]);
            }
            b[j + 1] = b[j + 1] - (divisor * b[pivotRow + 1]);
        }

        b[0] = b[0] - (functionDivisor * b[pivotRow + 1]);
        for(var p = 0; p < n; ++p) {
            c[p] = c[p] - (functionDivisor * A[pivotRow][p]);
        }
        b[pivotRow + 1] = b[pivotRow + 1] / pivot;
        for(var s = 0; s < n; ++s) {
            A[pivotRow][s] = A[pivotRow][s] / pivot;
        }
        printVector("C",c);
        printMatrix("A",A);
        printVector("B",b);
    }
    r=Result;
    r.A=A;
    r.c=c;
    r.b=b;
    r.pivotCol=pivotColumn;
    r.pivotRow=pivotRow;
    r.step=step;
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