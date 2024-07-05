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
    /**
     * ### 转向策略
     * 
     * - 左 "left"：一直会往左转。
     * - 右 "right"：一直会往右转。
     * - 黑 "black"：更倾向于链接黑色单元。
     * - 白 "white"：更倾向于链接白色单元。
     * - 少数派 "minority"：更倾向于链接在当前点给定范围内出现的最少的颜色。
     * - 多数派 "majority"：更倾向于链接出现的更多的颜色。
     * - 随机：随机转向。
     * 
     * 默认规则是少数派。 
     */
    turnPolicy: string;
    /**杂点判定阈值 */
    turdSize: number;
    optCurve: boolean;
    alphaMax: number;
    optTolerance: number;
    scale: number;
    /**填充还是描边 */
    isStroke: boolean;
}

export const DefaultConfig: Config = {
    isReady: false,
    turnPolicy: "minority",
    turdSize: 2,
    optCurve: true,
    alphaMax: 1,
    optTolerance: 0.2,
    scale: 1,
    isStroke: false,
};
