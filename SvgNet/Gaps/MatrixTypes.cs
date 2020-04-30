using System;
namespace SVG.Gaps {

    [Flags]
    internal enum MatrixTypes
    {
        Identity=0,
        Translation=1,
        Scaling=2,
        Unknown=4
    }
}
