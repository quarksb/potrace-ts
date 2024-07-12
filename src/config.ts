
/**
 * ### 参数配置
 * 
 * - `turnPolicy` 转向策略
 * - `turdSize` 杂点判定阈值
 * - `optCurve` 是否优化曲线
 * - `alphaMax` 角度阈值
 * - `optTolerance` 优化容差
 * - `scale` 缩放比例
 * - `isStroke` 填充还是描边
 */
export interface Config {
    /**
     * ### 转向策略
     * |Policy|Description|
     * |:--|:--|
     * |left:| which always takes a left turn |
     * |right:| which always takes a right turn |
     * |black:| which prefers to connect black components |
     * |white:| which prefers to connect white components |
     * |minority:| which prefers to connect the color (black or white) that occurs least frequently within a given neighborhood of the current position |
     * |majority:| which prefers to connect the color that occurs most frequently |
     * |random:| which makes a (more or less) random choice. |
     * 
     * The default turn policy is minority
     * 
     * @default "minority"
     */
    turnPolicy: string;
    /** turd 大便 turd size 杂点判定阈值, 面积小于此阈值的的路径将会被忽略 
     * @default 2
    */
    turdSize: number;
    /**是否优化曲线 
     * @default true
    */
    optCurve: boolean;
    /**弯曲度阈值, 用来判定是否是拐角 
     * @default 1
    */
    alphaMax: number;
    /**优化容差 
     * @default 0.2
    */
    optTolerance: number;
    /**缩放比例
     * @default 1
     */
    scale: number;
    /**填充还是描边 
     * @default false
    */
    isStroke: boolean;
}

export const DefaultConfig: Config = {
    turnPolicy: "minority",
    turdSize: 2,
    optCurve: true,
    alphaMax: 1,
    optTolerance: 0.2,
    scale: 1,
    isStroke: false,
};
