using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace MicroHomo
{
    public class MicroHomoInfo : GH_AssemblyInfo
    {
        public override string Name => "MicroHomo";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("507be107-2cdd-4925-b62f-a1d33a99403e");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";
    }
}