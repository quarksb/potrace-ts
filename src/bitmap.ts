import { Point } from "./base";

/**
 * Represents a bitmap image.
 */
export class Bitmap {
    /**
     * The width of the bitmap.
     */
    w: number;
    /**
     * The height of the bitmap.
     */
    h: number;
    /**
     * The total size of the bitmap (width * height).
     */
    size: number;
    /**
     * The Int8Array that represents the bitmap data.
     */
    data: Int8Array;

    /**
     * Creates a new Bitmap instance.
     * @param w The width of the bitmap.
     * @param h The height of the bitmap.
     */
    constructor (w: number, h: number) {
        this.w = w;
        this.h = h;
        this.size = w * h;
        this.data = new Int8Array(this.size);
    }

    /**
     * Checks if the specified coordinates are within the bitmap boundaries and the pixel at that position is set to 1.
     * @param x The x-coordinate.
     * @param y The y-coordinate.
     * @returns True if the pixel at the specified coordinates is set to 1, false otherwise.
     */
    checkPixel(x: number, y: number): boolean {
        return (
            x >= 0 && x < this.w &&
            y >= 0 && y < this.h &&
            this.data[this.w * y + x] === 1
        );
    }

    /**
     * Converts a linear index to a Point object representing the corresponding coordinates in the bitmap.
     * @param i The linear index.
     * @returns A Point object representing the coordinates in the bitmap.
     */
    index(i: number): Point {
        const y = Math.floor(i / this.w);
        const x = i % this.w;
        return new Point(x, y);
    }


    flip(x: number, y: number): void {
        this.data[this.w * y + x] = Number(!this.checkPixel(x, y));
    }

    /**
     * Creates a copy of the bitmap.
     * @returns A new Bitmap instance with the same dimensions and pixel values as the original bitmap.
     */
    copy(): Bitmap {
        const bm = new Bitmap(this.w, this.h);
        for (let i = 0; i < this.size; i++) {
            bm.data[i] = this.data[i];
        }
        return bm;
    }
}