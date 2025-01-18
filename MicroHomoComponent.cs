using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MN = MathNet.Numerics.LinearAlgebra.Double;
using MNL = MathNet.Numerics.LinearAlgebra;
using System.Linq;
using System.Diagnostics;

namespace MicroHomo
{
    public class MicroHomoComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public MicroHomoComponent()
          : base("MicroHomoComponent", "Homo",
            "Homogenize",
            "My", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMatrixParameter("dense", "dense", "dense", GH_ParamAccess.item);
            pManager.AddNumberParameter("E-mu-lx-ly", "E-mu-lx-ly", "E-mu-lx-ly", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMatrixParameter("Result", "Result", "Result", GH_ParamAccess.item);
            pManager.AddNumberParameter("Vol-Bulk-Shear-Possion", "Vol-Bulk-Shear-Possion", "Vol-Bulk-Shear-Possion", GH_ParamAccess.list);
            pManager.AddTextParameter("Timing", "Timing", "Timing", GH_ParamAccess.list);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Stopwatch totalStopwatch = Stopwatch.StartNew();

            Matrix matrixA = new Matrix(0, 0);
            List<double> ps = new List<double>();

            bool res1 = DA.GetData(0, ref matrixA);
            bool res2 = DA.GetDataList(1, ps);

            if (!res1 || !res2) { return; }


            double E = ps[0];
            double nu = ps[1];
            double lx = ps[2];
            double ly = ps[3];

            double tmp  = (E * nu) / ((1 + nu) * (1 - 2 * nu));
            double mu1 = E / (2 * (1 + nu));
            double lamb1 = (2 * tmp * mu1) / (tmp + 2 * mu1);
            double lamb2 = lamb1 * 1e-9;
            double mu2 = mu1 * 1e-9;


            int nelx = matrixA.ColumnCount;
            int nely = matrixA.RowCount;
            //int nelx = 100;
            //int nely = 100;
            double dx = lx / nelx;
            double dy = ly / nelx;
            int nel = nelx * nely;
            double a = 0.5 * dx;
            double b = 0.5 * dy;
            //double a = 0.005;
            //double b = 0.005;
            var CLambda = MN.DenseMatrix.OfArray(new double[,] { { 1, 1, 0 }, { 1, 1, 0 }, { 0, 0, 0 } });
            var CMu = MN.DenseMatrix .OfArray(new double[,] { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 1 } });
            var L = MN.DenseMatrix.OfArray(new double[,] { { 1, 0, 0 ,0 }, { 0, 0, 0, 1 }, { 0, 1, 1, 0 } });
            var keLambda = MN.DenseMatrix.Create(8, 8, 0.0);
            var keMu = MN.DenseMatrix.Create(8, 8, 0.0);
            var feLambda = MN.DenseMatrix.Create(8, 3, 0.0);
            var feMu = MN.DenseMatrix.Create(8, 3, 0.0);

            var dense = MN.DenseMatrix.Create(nely, nelx, 0.0);

            for (int i = 0; i < nely; i++)
            {
                for (int j = 0; j < nelx; j++)
                {
                    dense[i, j] = matrixA[i, j];
                }
            }
                //var dense = MN.DenseMatrix.OfArray(new double[,]{ { 0,1}, { 1, 0} });
                //var dense = MN.DenseMatrix.Create(nely, nelx, 1.0);
                //dense[0, 0] = 0.0;
                //dense[1, 1] = 0.0;
                //dense[2, 2] = 0.0;

            double[] xx = new double[2] { -0.577350269189626, 0.577350269189626 };
            double[] yy = xx;
            double[] ww = new double[2] { 1.0, 1.0 };

            Stopwatch Stopwatch1 = Stopwatch.StartNew();

            for (int ii = 0; ii < xx.Length; ++ii) 
            {  
                for (int jj =0; jj < yy.Length; ++jj)
                {
                    // Gauss points
                    double x = xx[ii];
                    double y = yy[jj];

                    // partial shape functions
                    double[] dNx = { 0.25 * (y - 1), 0.25 * (1 - y), 0.25 * (1 + y), -0.25 * (1 + y) };
                    double[] dNy = { 0.25 * (x - 1), -0.25 * (1 + x), 0.25 * (1 + x), 0.25 * (1 - x) };

                    var dNtmp = MN.DenseMatrix.OfArray(new double[,] { { 0.25 * (y - 1), 0.25 * (1 - y), 0.25 * (1 + y), -0.25 * (1 + y) },
                                                                       { 0.25 * (x - 1), -0.25 * (1 + x), 0.25 * (1 + x), 0.25 * (1 - x) } });

                    var jtmp = MN.DenseMatrix.OfArray(new double[,] { { -a, -b }, { a, -b}, { a, b }, {- a , b} });
                    var J = dNtmp * jtmp;
                    double detJ  = J.Determinant();
                    var invJ = J.Inverse();
                    double weight = ww[ii] * ww[jj] * detJ;

                    var G = MN.DenseMatrix.Create(4, 4, 0.0);
                    G.SetSubMatrix(0, 0, invJ);
                    G.SetSubMatrix(2, 2, invJ);

                    var dN = MN.DenseMatrix.Create(4, 8, 0.0);
                    dN[0, 0] = dNx[0]; dN[0, 2] = dNx[1]; dN[0, 4] = dNx[2]; dN[0, 6] = dNx[3];
                    dN[1, 0] = dNy[0]; dN[1, 2] = dNy[1]; dN[1, 4] = dNy[2]; dN[1, 6] = dNy[3];
                    dN[2, 1] = dNx[0]; dN[2, 3] = dNx[1]; dN[2, 5] = dNx[2]; dN[2, 7] = dNx[3];
                    dN[3, 1] = dNy[0]; dN[3, 3] = dNy[1]; dN[3, 5] = dNy[2]; dN[3, 7] = dNy[3];

                    var B = L * G * dN;

                    var e0 = MN.DenseMatrix.CreateIdentity(3);
                    keLambda.Add(weight * B.Transpose() * CLambda * B, keLambda);
                    keMu.Add(weight * B.Transpose() * CMu * B, keMu);
                    feLambda.Add(weight * B.Transpose() * CLambda * e0, feLambda);
                    feMu.Add(weight * B.Transpose() * CMu * e0, feMu);
                }
            }

            Stopwatch1.Stop();
            string Time1 = $"Element stiffness matrix assamble Time: {Stopwatch1.ElapsedMilliseconds} ms";

            Stopwatch Stopwatch2 = Stopwatch.StartNew();

            var nodenrs = MN.DenseMatrix.OfColumnMajor(nely + 1, nelx + 1, Enumerable.Range(0, (nely + 1) * (nelx + 1)).Select(x => (double)x));
            double[] edofVec = (nodenrs.SubMatrix(0, nely, 0, nelx) * 2 + 2).AsColumnMajorArray();
            //int[,] nodenrs = new int[1 + nely, 1 + nelx];
            //int count = 0;
            //for (int i = 0; i <= nely; i++)
            //{
            //    for (int j = 0; j <= nelx; j++)
            //    {
            //        nodenrs[j, i] = count++;
            //    }
            //}

            //int[] edofVec = new int[nel];
            //count = 0;
            //for (int i = 0; i < nely; i++)
            //{
            //    for (int j = 0; j < nelx; j++)
            //    {
            //        edofVec[count++] = 2 * nodenrs[j, i] + 2;
            //    }
            //}

            var edofMat = MN.DenseMatrix.OfArray(new double[nel, 8]);
            var indexTmp = MN.DenseVector.OfArray(new double[] { 0, 1, 2 * nely + 2, 2 * nely + 3, 2 * nely, 2 * nely + 1, -2, -1 });
            for (int i =0; i<nel; i++)
            {
                edofMat.SetRow(i, edofVec[i] + indexTmp);
            }

            //for (int i = 0; i < nel; i++)
            //{
            //    for (int j = 0; j < 8; j++)
            //    {
            //        edofMat[i, j] = edofVec[i] + new[] { 0, 1, 2 * nely + 2, 2 * nely + 3, 2 * nely, 2 * nely + 1, -2, -1 }[j];
            //    }
            //}

            int nn = (nelx + 1) * (nely + 1);
            int nnP = nelx * nely;
            var nnPArray = MN.DenseMatrix.OfArray(new double[nely + 1, nelx + 1]);
            //for (int i = 0; i < nelx; i++)
            //{
            //    for (int j = 0; j < nely; j++)
            //    {
            //        nnPArray[j, i] = i * nelx + j ;
            //    }
            //}

            var nnPArray0 = MN.DenseMatrix.OfColumnMajor(nely, nelx, Enumerable.Range(0, nnP).Select(x => (double)x).ToArray());
            nnPArray.SetSubMatrix(0, 0, nnPArray0);
            nnPArray.SetRow(nely, nnPArray.Row(0));
            nnPArray.SetColumn(nelx, nnPArray.Column(0));
            //for (int j = 0; j < nelx; j++)
            //{
            //    nnPArray[nely, j] = nnPArray[0, j];
            //}

            //for (int i = 0; i < nely; i++)
            //{
            //    nnPArray[i, nelx] = nnPArray[i, 0];
            //}

            var dofVector = new int[2 * nn];
            for (int i = 0; i < nelx + 1; i++)
            {
                for (int j = 0; j < nely + 1; j++)
                {
                    dofVector[2 * (j + (nely + 1) * i)] = 2 * (int)nnPArray[j, i];
                    dofVector[2 * (j + (nely + 1) * i) + 1] = 2 * (int)nnPArray[j, i] + 1;
                }
            }

            for (int i = 0; i < nel; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    edofMat[i, j] = dofVector[(int)edofMat[i, j]];
                }
            }

            int ndof = 2 * nnP;

            // Assembly global stiffness matrix
            //var iK_ = new MN.DenseMatrix(nel * 8, 8);
            //var jK_ = new MN.DenseMatrix(nel, 8 * 8);
            //for (int i = 0; i < nel; i++)
            //{
            //    for (int j = 0; j < 8; j++)
            //    {
            //        for (int k = 0; k < 8; k++)
            //        {
            //            iK_[i * 8 + k, j] = edofMat[i, j];
            //            jK_[i, j * 8 + k] = edofMat[i, j];
            //        }
            //    }
            //}

            var iK = edofMat.KroneckerProduct(MN.DenseMatrix.Create(8, 1, 1.0)).Transpose();
            var jK = edofMat.KroneckerProduct(MN.DenseMatrix.Create(1, 8, 1.0)).Transpose();
            //var iK = iK_.Transpose();
            //var jK = jK_.Transpose();

            var lambda = new MN.DenseMatrix(nely, nelx);
            var mu = new MN.DenseMatrix(nely, nelx);
            for (int i = 0; i < nely; i++)
            {
                for (int j = 0; j < nelx; j++)
                {
                    //if (dense[i, j] > 0.5)
                    //{
                    //    lambda[i, j] = lamb1;
                    //    mu[i, j] = mu1;
                    //}
                    //else
                    //{
                    //    lambda[i, j] = lamb2;
                    //    mu[i, j] = mu2;
                    //}
                    lambda[i, j] = dense[i, j] > 0.5 ? lamb1 : lamb2;
                    mu[i, j] = dense[i, j] > 0.5 ? mu1 : mu2;
                }
            }

            var keLambdaVec = MN.DenseVector.OfArray(keLambda.ToColumnMajorArray());
            var keMuVec = MN.DenseVector.OfArray(keMu.ToColumnMajorArray());
            var lambdaVec = MN.DenseVector.OfArray(lambda.ToColumnMajorArray());
            var muVec = MN.DenseVector.OfArray(mu.ToColumnMajorArray());

            //var sK = new MN.DenseMatrix(64, nel);

            var sK = keLambdaVec.OuterProduct(lambdaVec) + keMuVec.OuterProduct(muVec);
            //for (int i = 0; i < keLambdaVec.Count; i++)
            //{
            //    for (int j = 0; j < lambdaVec.Count; j++)
            //    {
            //        sK[i, j] = keLambdaVec[i] * lambdaVec[j] + keMuVec[i] * muVec[j];
            //    }
            //}

            var iKarray = iK.ToColumnMajorArray();
            var jKarray = jK.ToColumnMajorArray();
            var sKarray = sK.ToColumnMajorArray();

            //int[] iKidx = new int[iKarray.Length];
            //int[] jKidx = new int[jKarray.Length];
            //for (int i = 0; i < iKarray.Length; i++)
            //{
            //    iKidx[i] = (int)iKarray[i];
            //    jKidx[i] = (int)jKarray[i];
            //}
            List<int> iKidx = new List<int>();
            List<int> jKidx = new List<int>();
            List<double> sKval = new List<double>();
            for (int i = 0; i < iKarray.Length; i++)
            {
                if (!((int)iKarray[i] == 0)|| !((int)iKarray[i] == 1) || !((int)jKarray[i] == 0) || !((int)jKarray[i] == 1))
                {
                    iKidx.Add((int)iKarray[i]);
                    jKidx.Add((int)jKarray[i]);
                    sKval.Add(sKarray[i]);
                }
            }
            iKidx.Add(0);
            iKidx.Add(1);
            jKidx.Add(0);
            jKidx.Add(1);
            sKval.Add(1.0);
            sKval.Add(1.0);



            //List<Tuple<int, int, double>> vals = new List<Tuple<int, int, double>>();

            //for (int i = 0; i < iKarray.Length; i++)
            //{
            //    vals.Add(new Tuple<int, int, double>((int)iKarray[i], (int)jKarray[i], sKarray[i]));
            //}

            //var K = MN.SparseMatrix.OfIndexed(ndof, ndof, vals);

            var K = MNL.Matrix<double>.Build.SparseFromCoordinateFormat(ndof, ndof, sKarray.Length, iKidx.ToArray(), jKidx.ToArray(), sKval.ToArray());

            Stopwatch2.Stop();
            string Time2 = $"Global stiffness matrix assamble Time: {Stopwatch2.ElapsedMilliseconds} ms";

            Stopwatch Stopwatch3 = Stopwatch.StartNew();

            var feLambdaVec = MN.DenseVector.OfArray(feLambda.ToColumnMajorArray());
            var feMuVec = MN.DenseVector.OfArray(feMu.ToColumnMajorArray());
            //var sF = new MN.DenseMatrix(24, nel);
            var iF = new MN.DenseMatrix(8 * 3, nel);
            var jF = new MN.DenseMatrix(8 * 3, nel);

            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < nel; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        iF[i + 8 * k, j] = edofMat[j, i];
                    }
                }
            }

            jF.SetSubMatrix(0, 0, MN.DenseMatrix.Create(8, nel, 0));
            jF.SetSubMatrix(1 * 8, 0, MN.DenseMatrix.Create(8, nel, 1));
            jF.SetSubMatrix(2 * 8, 0, MN.DenseMatrix.Create(8, nel, 2));

            //for (int i = 0; i < feLambdaVec.Count; i++)
            //{
            //    for (int j = 0; j < lambdaVec.Count; j++)
            //    {
            //        sF[i, j] = feLambdaVec[i] * lambdaVec[j] + feMuVec[i] * muVec[j];
            //    }
            //}
            var sF = feLambdaVec.OuterProduct(lambdaVec) + feMuVec.OuterProduct(muVec);

            var iFarray = iF.ToColumnMajorArray();
            var jFarray = jF.ToColumnMajorArray();
            var sFarray = sF.ToColumnMajorArray();

            int[] iFidx = new int[iFarray.Length];
            int[] jFidx = new int[jFarray.Length];

            for (int i = 0; i < iFarray.Length; i++)
            {
                iFidx[i] = (int)iFarray[i];
                jFidx[i] = (int)jFarray[i];
            }

            var F = MNL.Matrix<double>.Build.SparseFromCoordinateFormat(ndof, 3, sFarray.Length, iFidx, jFidx, sFarray);
            F.SetRow(0, new double[] { 0.0, 0.0, 0.0 });
            F.SetRow(1, new double[] { 0.0, 0.0, 0.0 });

            Stopwatch3.Stop();
            string Time3 = $"External force assamble Time: {Stopwatch3.ElapsedMilliseconds} ms";

            Stopwatch Stopwatch4 = Stopwatch.StartNew();

            var solver = new MN.Solvers.BiCgStab();

            double[] zeroTmp = Enumerable.Repeat(0.0, ndof).ToArray();
            //K.SetRow(0, zeroTmp);
            //K.SetRow(1, zeroTmp);
            //K.SetColumn(0, zeroTmp);
            //K.SetColumn(1, zeroTmp);
            //K[0, 0] = 1.0;
            //K[1, 1] = 1.0;

            //var Farray = F.AsColumnArrays();
            //var F1 = MN.SparseVector.OfVector(F.Column(0));
            //var F2 = MN.SparseVector.OfVector(F.Column(1));
            //var F3 = MN.SparseVector.OfVector(F.Column(2));
            //F.SetRow(0, new double[] { 1, 1, 1 });

            //var chiLoc0 = K.SolveIterative(MN.SparseVector.OfVector(F.Column(0)), solver);
            //var chiLoc1 = K.SolveIterative(MN.SparseVector.OfVector(F.Column(1)), solver);
            //var chiLoc2 = K.SolveIterative(MN.SparseVector.OfVector(F.Column(2)), solver);
            //var chi = MN.DenseMatrix.Create(ndof, 3, 0.0);
            //chi.SetColumn(0, chiLoc0);
            //chi.SetColumn(1, chiLoc1);
            //chi.SetColumn(2, chiLoc2);
            var preconditioner = new MN.Solvers.DiagonalPreconditioner();
            var chi = K.SolveIterative(F, solver, preconditioner);
            Stopwatch4.Stop();
            string Time4 = $"Solve Time: {Stopwatch4.ElapsedMilliseconds} ms";

            ////var chiLoc0 = K.SubMatrix(2, ndof - 2, 2, ndof - 2).SolveIterative(MN.DenseVector.OfArray(Farray[0]), solver);
            ////var chiLoc1 = K.SubMatrix(2, ndof - 2, 2, ndof - 2).SolveIterative(MN.DenseVector.OfArray(Farray[1]), solver);
            ////var chiLoc2 = K.SubMatrix(2, ndof - 2, 2, ndof - 2).SolveIterative(MN.DenseVector.OfArray(Farray[2]), solver);
            //////var ttt = F.SubMatrix(2, ndof - 2, 0, 1).ToColumnMajorArray();
            //var chi = MN.DenseMatrix.Create(ndof, 3, 0.0);
            //chi.SetColumn(0, chiLoc0);
            //chi.SetColumn(1, chiLoc1);
            //chi.SetColumn(2, chiLoc2);

            Stopwatch Stopwatch5 = Stopwatch.StartNew();

            var ke = keMu + keLambda;
            var fe = feMu + feLambda;
            zeroTmp = Enumerable.Repeat(0.0, 8).ToArray();
            ke.SetRow(0, zeroTmp);
            ke.SetRow(1, zeroTmp);
            ke.SetRow(3, zeroTmp);
            ke.SetColumn(0, zeroTmp);
            ke.SetColumn(1, zeroTmp);
            ke.SetColumn(3, zeroTmp);
            ke[0, 0] = 1.0;
            ke[1, 1] = 1.0;
            ke[3, 3] = 1.0;
            fe.SetRow(0, new double[] { 0.0, 0.0, 0.0});
            fe.SetRow(1, new double[] { 0.0, 0.0, 0.0 });
            fe.SetRow(3, new double[] { 0.0, 0.0, 0.0 });

            var chi0_e = ke.SolveIterative(fe, solver, preconditioner);
            //var chi0_e1 = ke.SolveIterative(fe.Column(1), solver);
            //var chi0_e2 = ke.SolveIterative(fe.Column(2), solver);

            //var chi0_0 = MN.DenseMatrix.Create(nel, 8, 0.0);
            //var chi0_1 = MN.DenseMatrix.Create(nel, 8, 0.0);
            //var chi0_2 = MN.DenseMatrix.Create(nel, 8, 0.0);

            //for (int i = 0; i < nel; i++)
            //{
            //    chi0_0.SetRow(i, chi0_e.Column(0));
            //    chi0_1.SetRow(i, chi0_e.Column(1));
            //    chi0_2.SetRow(i, chi0_e.Column(2));
            //}

            var chi0_0 = chi0_e.Column(0).OuterProduct(MN.DenseVector.Create(nel, 1.0)).Transpose();
            var chi0_1 = chi0_e.Column(1).OuterProduct(MN.DenseVector.Create(nel, 1.0)).Transpose();
            var chi0_2 = chi0_e.Column(2).OuterProduct(MN.DenseVector.Create(nel, 1.0)).Transpose();
            List<MNL.Matrix<double>> chi0 = new List<MNL.Matrix<double>> { chi0_0, chi0_1, chi0_2 };

            var CH = MN.DenseMatrix.Create(3, 3, 0.0);
            double cellVolume = lx * ly;

            for (int i = 0; i< 3; i++)
            {
                for (int j = 0; j < 3; j++) 
                {
                    var chiTmpi = MN.DenseMatrix.Create(nel, 8, 0.0);
                    var chiTmpj = MN.DenseMatrix.Create(nel, 8, 0.0);

                    for (int ii = 0; ii < nel; ii++)
                    {
                        for (int jj = 0; jj < 8; jj++) 
                        {
                            chiTmpi[ii, jj] = chi.AsColumnMajorArray()[(int)edofMat[ii, jj] + (i * ndof)];
                        }
                    }

                    for (int ii = 0; ii < nel; ii++)
                    {
                        for (int jj = 0; jj < 8; jj++)
                        {
                            chiTmpj[ii, jj] = chi.AsColumnMajorArray()[(int)edofMat[ii, jj] + (j * ndof)];
                        }
                    }

                    var leftTmp = chi0[i] - chiTmpi;
                    var rightTmp = chi0[j] - chiTmpj;
                    var sumLambda_ = (leftTmp * keLambda).PointwiseMultiply(rightTmp);
                    var sumMu_ = (leftTmp * keMu).PointwiseMultiply(rightTmp);
                    var sumLambdaArray = sumLambda_.RowSums().ToArray();
                    var sumMuArray = sumMu_.RowSums().ToArray();
                    var sumLambda = MN.DenseMatrix.OfColumnMajor(nely, nelx, sumLambdaArray);
                    var sumMu = MN.DenseMatrix.OfColumnMajor(nely, nelx, sumMuArray);

                    var onesTmp = MN.DenseVector.Create(nelx, 1.0);

                    CH[i, j] = (1 / cellVolume) * (onesTmp * (sumLambda.PointwiseMultiply(lambda) + sumMu.PointwiseMultiply(mu)) * onesTmp);
                }
            }
            Stopwatch5.Stop();
            string Time5 = $"Other Time: {Stopwatch5.ElapsedMilliseconds} ms";

            totalStopwatch.Stop();
            string totalTime = $"Total Time: {totalStopwatch.ElapsedMilliseconds} ms";

            double vol = dense.Enumerate().Sum() / (nelx * nely); ;

            Matrix result = new Matrix(3, 3);

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    result[i, j] = CH[i, j];
                }
            }

            double bulkModulus = (CH[0, 0] + 2 * CH[0, 1]) * 0.3333333333333333333333;
            double shearModulus = CH[2, 2];
            double possionRatio = CH[0, 1] / (CH[0, 0] + CH[0, 1]);

            DA.SetData(0, result);
            DA.SetDataList(1,  new List<double> { vol, bulkModulus, shearModulus, possionRatio } );
            DA.SetDataList(2, new List<string> { Time1, Time2, Time3, Time4, Time5, totalTime } );
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// You can add image files to your project resources and access them like this:
        /// return Resources.IconForThisComponent;
        /// </summary>
        protected override System.Drawing.Bitmap Icon => null;

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid => new Guid("5c044db3-e148-4f65-84e5-18e73ebb4be3");
    }
}