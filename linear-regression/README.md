## Linear Regression

### Usage
Perform linear regression of the form:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y = x*b + e

By solving a linear programming problem of the form:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; min&nbsp;&nbsp;&nbsp;&nbsp;|y - x*b|

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; s.t.&nbsp;&nbsp;&nbsp;&nbsp; -M <= b <= M

### Variables

    solution = linearRegression(x,y)
    solution = linearRegression(x,y,M,opt_params)

#### Inputs

*n* is the number of samples and *m* is the number of explanatory variables 

| Variable Name | Dimension      | Description                                                                                                              | Required? |
| ------------- | -------------- | ---------------------                                                                                                    | --------- |
| x             | n x m          | Explanatory variables                                                                                                    | Yes       |
| y             | n x 1 or n x m | Response variable(s)                                                                                                     | Yes       |
| M             | m+1 x m        |                                                                                                                          | No        |
| opt_params    | n/a            | Structure containing [Gurobi parameters](https://www.gurobi.com/documentation/8.1/refman/parameters.html#sec:Parameters) | No        |

#### Output

solution is a structure

| Field Name | Dimension          | Description              |
| ---------- | ------------------ | ------------------------ |
| B          | m+1 x 1 or m+1 x m | Regression coefficients  |
| w          | n x m              | Loss                     |
| y          | n x 1 or n x m     | Estimated y values (x*b) |
| status     | n/a                | Optimization status      |

## Subset Selection

### Usage
Perform linear regression of the form:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y = x*b + e

By solving a linear programming problem of the form:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; min&nbsp;&nbsp;&nbsp;&nbsp;|y - x*b|

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; s.t.&nbsp;&nbsp;&nbsp;&nbsp; -Mz <= b <= Mz
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; sum(z) <= p

### Variables

    solution = bssp(x,y)
    solution = bssp(x,y,p,M,opt_params)

#### Inputs

*n* is the number of samples and *m* is the number of explanatory variables 

| Variable Name | Dimension      | Description                                                                                                              | Required? |
| ------------- | -------------- | ---------------------                                                                                                    | --------- |
| x             | n x m          | Explanatory variables                                                                                                    | Yes       |
| y             | n x 1 or n x m | Response variable(s)                                                                                                     | Yes       |
| p             | scalar       | Number of predictors                                                                                                     | No        |
| M             | m+1 x m        |                                                                                                                          | No        |
| opt_params    | n/a            | Structure containing [Gurobi parameters](https://www.gurobi.com/documentation/8.1/refman/parameters.html#sec:Parameters) | No        |

#### Output

solution is a structure

| Field Name | Dimension          | Description                                     |
| ---------- | ------------------ | ----------------------------------------------- |
| B          | m+1 x 1 or m+1 x m | Regression coefficients                         |
| w          | n x m              | Loss                                            |
| z          | m+1 x 1            | Indicates if predictor is active (1) or not (0) |
| y          | n x 1 or n x m     | Estimated y values (x*b)                        |
| status     | n/a                | Optimization status                             |

## Dependencies
Requires the [Gurobi](https://www.gurobi.com/) optimizer

## License
[MIT](https://choosealicense.com/licenses/mit/)
