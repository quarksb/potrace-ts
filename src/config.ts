/*parameters:
 *        turnpolicy ("black" / "white" / "left" / "right" / "minority" / "majority")
 *          how to resolve ambiguities in path decomposition. (default: "minority")       
 *        turdsize
 *          suppress speckles of up to this size (default: 2)
 *        optcurve (true / false)
 *          turn on/off curve optimization (default: true)
 *        alphamax
 *          corner threshold parameter (default: 1)
 *        opttolerance 
 *          curve optimization tolerance (default: 0.2)
 */
export interface Config {
    isReady: boolean;
    /**转向策略 */
    turnPolicy: string;
    /**杂点判定阈值 */
    turdSize: number;
    optCurve: boolean;
    alphaMax: number;
    optTolerance: number;
}

export const DefaultConfig: Config = {
    isReady: false,
    turnPolicy: "minority",
    turdSize: 2,
    optCurve: true,
    alphaMax: 1,
    optTolerance: 0.2
};
