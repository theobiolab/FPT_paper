(* ::Package:: *)

Print[$CommandLine]

f = $CommandLine[[4]];
customCondition = ToExpression[$CommandLine[[5]]];
conditionString = StringReplace[$CommandLine[[5]], {" " -> ""}];
outfname = "./results/results_lt0_" <> First[StringSplit[$CommandLine[[4]], "_"]] <> "_" <> conditionString <> ".csv"; (* Construct output CSV file name *)
timec = 600;
outstream = OpenWrite[outfname];
WriteString[outstream, "Line, Expression, Evaluation\n"]; (* Write CSV header *)
file = OpenRead[f];
line = ReadLine[file];
i = 1;
While[Not[SameQ[line, EndOfFile]],
    coefexpr = StringSplit[line, ";"];
    expr = ToExpression[coefexpr[[2]]];
    (*
    out1 = With[{assum = customCondition}, 
                TimeConstrained[Simplify[Reduce[assum \[Implies] expr < 0, Reals], assum], timec]];
    *)
    out1 = Assuming[customCondition, 
         TimeConstrained[
           Simplify[Reduce[expr < 0, Reals], customCondition], 
           timec
         ]
       ];  

    If[out1 === True || out1 === False,
        evaluation = If[out1 === True, "True", "False"],

        out1 = With[{assum = customCondition}, Simplify[Reduce[assum \[Implies] out1, Reals], assum]]
    ]
    
    WriteString[outstream, ToString[i] <> ", " <> ToString[$CommandLine[[5]]] <> ", " <> coefexpr[[2]] <> "<0, " <> ToString[out1] <> "\n"];
    line = ReadLine[file];
    i = i + 1;
];
Close[file];
Close[outstream];

