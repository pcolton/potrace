//
//  Potrace.swift
//  PhotoHero
//
//  Created by Paul Colton on 5/19/17.
//  Copyright © 2017 Paul Colton. All rights reserved.
//
//  Licensed under the GPL v2.0
//
//  See https://github.com/pcolton/potrace
//
//  Swift port of https://github.com/kilobtye/potrace/blob/master/potrace.js
//

import Foundation
import UIKit

final class Potrace {

    class Point<T> {
        var x: T
        var y: T
        
        init(x: T, y: T) {
            self.x = x
            self.y = y
        }
        
        func copy() -> Point<T> {
            return Point<T>(x: self.x, y: self.y)
        }
        
        static func toPointF(point: PointI) -> PointF {
            return PointF(x: Double(point.x), y: Double(point.y))
        }
    }
    
    typealias PointF = Point<Double>
    typealias PointI = Point<Int>
    
    struct Bitmap {
        var w: Int
        var h: Int
        var size: Int
        var data: [UInt8]

        init(width: Int, height: Int) {
            self.w = width
            self.h = height
            self.size = width * height
            self.data = [UInt8](repeating: 0, count: self.size)
        }
        
        func at(x: Int, y: Int) -> Bool {
            return (x >= 0 && x < self.w && y >= 0 && y < self.h) &&
                data[self.w * y + x] == 1
        }
        
        func index(i: Int) -> PointI {
            let point = PointI(x: 0, y: 0)
            point.y = Int(floor(Double(i / self.w)))
            point.x = i - point.y * self.w
            return point
        }
        
        mutating func flip(x: Int, y: Int) {
            if self.at(x: x, y: y) {
                data[self.w * y + x] = 0
            } else {
                data[self.w * y + x] = 1
            }
        }
        
        func copy() -> Bitmap {
            var bm = Bitmap(width: w, height: h)
            for (index, value) in data.enumerated() {
                bm.data[index] = value
            }
            return bm
        }
    }

    class Curve {
        var n: Int
        var tag: [String]
        var c: BoxedArray<PointF>
        var alphacurve: Int
        var vertex: BoxedArray<PointF>
        var alpha: [Double]
        var alpha0: [Double]
        var beta: [Double]
        
        init(n: Int) {
            self.n = n
            self.tag = [String](repeating: "", count: n)
            self.c = BoxedArray(repeating: Point(x: 0.0, y: 0.0), count: n * 3)
            self.alphacurve = 0
            self.vertex = BoxedArray(repeating: Point(x: 0.0, y: 0.0), count: n)
            self.alpha = [Double](repeating: 0.0, count: n)
            self.alpha0 = [Double](repeating: 0.0, count: n)
            self.beta = [Double](repeating: 0.0, count: n)
        }
    }
    
    struct Sum {
        var x: Int
        var y: Int
        var xy: Int
        var x2: Int
        var y2: Int
    }

    class Path {
        var area: Int
        var len: Int
        var curve: Curve
        var pt: BoxedArray<PointI>
        var minX: Int
        var minY: Int
        var maxX: Int
        var maxY: Int
        var sign: String
        var sums: BoxedArray<Sum>
        var lon: [Int]
        var po: [Int]
        var m: Int
        var x0: Int
        var y0: Int
        
        init() {
            self.area = 0
            self.len = 0
            self.curve = Curve(n: 0)
            self.pt = BoxedArray()
            self.minX = 100000
            self.minY = 100000
            self.maxX = -1
            self.maxY = -1
            self.sign = ""
            self.sums = BoxedArray()
            self.lon = []
            self.po = []
            self.m = 0
            self.x0 = 0
            self.y0 = 0
        }
    }
    
    struct Info {
        var isReady: Bool
        var turnpolicy: String
        var turdsize: Int
        var optcurve: Bool
        var alphamax: Double
        var opttolerance: Double
        
        init() {
            self.isReady = false
            self.turnpolicy = "minority"
            self.turdsize = 2
            self.optcurve = true
            self.alphamax = 1.0
            self.opttolerance = 0.2
        }
        
    }
    
    fileprivate var bm: Bitmap!
    fileprivate var pathlist = [Path]()
    fileprivate var info = Info()
    
    init(data: UnsafeMutableRawPointer, width: Int, height: Int) {
        let pixelData = data.assumingMemoryBound(to: UInt8.self)

        bm = Bitmap(width: width, height: height)

        var j = 0, i = 0
        
        for _ in 0..<bm.data.count {
            let val = 0.2126 * Double(pixelData[i]) + 0.7153 * Double(pixelData[i + 1]) +
                0.0721 * Double(pixelData[i + 2])
            bm.data[j] = val < 128 ? 1 : 0
            i += 4
            j += 1
        }
    }

    func process() {
        bmToPathList()
        processPath()
    }
    
    func bmToPathList() {
        var bm1 = bm.copy()
        var currentPoint = Point(x: 0, y: 0)
        var path: Path
        
        func findNext(point: PointI) -> PointI? {
            var i = bm1.w * point.y + point.x
            while i < bm1.size && bm1.data[i] != 1 {
                i += 1
            }
            if i < bm1.size {
                return bm1.index(i: i)
            } else {
                return nil
            }
        }
        
        func majority(x: Int, y: Int) -> Int {

            for i in 2..<5 {
                var ct = 0

                for a in (-i + 1)...(i - 1) {
                    ct += bm1.at(x: x + a, y: y + i - 1) ? 1 : -1
                    ct += bm1.at(x: x + i - 1, y: y + a - 1) ? 1 : -1
                    ct += bm1.at(x: x + a - 1, y: y - i) ? 1 : -1
                    ct += bm1.at(x: x - i, y: y + a) ? 1 : -1
                }
                
                return ct > 0 ? 1 : 0
            }
            
            return 0
        }

        func findPath(point: PointI) -> Path {
            var path = Path(),
            x = point.x, y = point.y,
            dirx = 0, diry = 1, tmp: Int
            
            path.sign = bm.at(x: point.x, y: point.y) ? "+" : "-"
            
            while true {
                path.pt.append(Point(x: x, y: y))
                if x > path.maxX {
                    path.maxX = x
                }
                if x < path.minX {
                    path.minX = x
                }
                if y > path.maxY {
                    path.maxY = y
                }
                if y < path.minY {
                    path.minY = y
                }
                path.len += 1
                
                x += dirx
                y += diry
                path.area -= x * diry
                
                if x == point.x && y == point.y {
                    break
                }
                
                let l = bm1.at(x: x + (dirx + diry - 1 ) / 2, y: y + (diry - dirx - 1) / 2)
                let r = bm1.at(x: x + (dirx - diry - 1) / 2, y: y + (diry + dirx - 1) / 2)
                
                if r && !l {
                    if info.turnpolicy == "right" ||
                        (info.turnpolicy == "black" && path.sign == "+") ||
                        (info.turnpolicy == "white" && path.sign == "-") ||
                        (info.turnpolicy == "majority" && majority(x: x, y: y) != 0) ||
                        (info.turnpolicy == "minority" && majority(x: x, y: y) == 0) {
                        tmp = dirx
                        dirx = -diry
                        diry = tmp
                    } else {
                        tmp = dirx
                        dirx = diry
                        diry = -tmp
                    }
                } else if r {
                    tmp = dirx
                    dirx = -diry
                    diry = tmp
                } else if !l {
                    tmp = dirx
                    dirx = diry
                    diry = -tmp
                }
            }
            
            return path
        }

        func xorPath(path: Path) {
            var y1 = path.pt[0].y,
            len = path.len,
            x: Int, y: Int, maxX: Int, minY: Int
            
            for i in 1..<len {
                x = path.pt[i].x
                y = path.pt[i].y
                
                if y != y1 {
                    minY = y1 < y ? y1 : y
                    maxX = path.maxX
                    for j in x..<maxX {
                        bm1.flip(x: j, y: minY)
                    }
                    y1 = y
                }
            }
        }
        
        while let currentPoint = findNext(point: currentPoint) {
            
            path = findPath(point: currentPoint)
            
            xorPath(path: path)
            
            if path.area > info.turdsize {
                pathlist.append(path)
            }
        }
    }
    
    //
    // processPath
    //

    func processPath() {
        
        struct Quad {
            var data: [Double]
            
            init() {
                self.data = [Double](repeating: 0.0, count: 9)
            }
            
            func at(x: Int, y: Int) -> Double {
                return self.data[x * 3 + y]
            }
        }
        
        func mod(a: Int, n: Int) -> Int {
            return a >= n ? a % n : a>=0 ? a : n-1-(-1-a) % n
        }

        func xprod(p1: PointI, p2: PointI) -> Int {
            return p1.x * p2.y - p1.y * p2.x
        }
        
        func cyclic(a: Int, b: Int, c: Int) -> Bool {
            if a <= c {
                return (a <= b && b < c)
            } else {
                return (a <= b || b < c)
            }
        }
        
        func sign(i: Int) -> Int {
            return i > 0 ? 1 : i < 0 ? -1 : 0
        }
        
        func quadform(Q: Quad, w: PointF) -> Double {
            var v = [Double](repeating: 0.0, count: 3)
            var sum: Double = 0.0
            
            v[0] = w.x
            v[1] = w.y
            v[2] = 1
            
            for i in 0..<3 {
                for j in 0..<3 {
                    sum += Double(v[i]) * Q.at(x: i, y: j) * Double(v[j])
                }
            }
            return sum
        }
        
        func interval(lambda: Double, a: PointF, b: PointF) -> PointF {
            let res = PointF(x: 0, y: 0)
            
            res.x = a.x + lambda * (b.x - a.x)
            res.y = a.y + lambda * (b.y - a.y)
            return res
        }
        
        func dorth_infty(p0: PointF, p2: PointF) -> PointI {
            let r = PointI(x: 0, y: 0)
            
            r.y = sign(i: Int(p2.x - p0.x))
            r.x = -sign(i: Int(p2.y - p0.y))
            
            return r
        }
        
        func ddenom(p0: PointF, p2: PointF) -> Double {
            let r = dorth_infty(p0: p0, p2: p2)
            
            return Double(r.y) * (p2.x - p0.x) - Double(r.x) * (p2.y - p0.y)
        }
        
        func dpara(p0: PointF, p1: PointF, p2: PointF) -> Double {

            let x1 = p1.x - p0.x
            let y1 = p1.y - p0.y
            let x2 = p2.x - p0.x
            let y2 = p2.y - p0.y
            
            return x1 * y2 - x2 * y1
        }
        
        func cprod(p0: PointF, p1: PointF, p2: PointF, p3: PointF) -> Double {

            let x1 = p1.x - p0.x
            let y1 = p1.y - p0.y
            let x2 = p3.x - p2.x
            let y2 = p3.y - p2.y
            
            return x1 * y2 - x2 * y1
        }
        
        func iprod(p0: PointF, p1: PointF, p2: PointF) -> Double {

            let x1 = p1.x - p0.x
            let y1 = p1.y - p0.y
            let x2 = p2.x - p0.x
            let y2 = p2.y - p0.y
            
            return x1*x2 + y1*y2
        }
        
        func iprod1(p0: PointF, p1: PointF, p2: PointF, p3: PointF) -> Double {
            
            let x1 = p1.x - p0.x
            let y1 = p1.y - p0.y
            let x2 = p3.x - p2.x
            let y2 = p3.y - p2.y
            
            return x1 * x2 + y1 * y2
        }
        
        func ddist(p: PointF, q: PointF) -> Double {
            let part1 = (p.x - q.x) * (p.x - q.x)
            let part2 = (p.y - q.y) * (p.y - q.y)
            return sqrt(part1 + part2)
        }
        
        func bezier(t: Double, p0: PointF, p1: PointF, p2: PointF, p3: PointF) -> PointF {
            let s = 1.0 - t
            let res = PointF(x: 0, y: 0)
            let part1 = s*s*s*p0.x
            let part11 = 3*(s*s*t)*p1.x
            let part2 = 3*(t*t*s)*p2.x
            let part22 = t*t*t*p3.x
            let part3 = s*s*s*p0.y
            let part33 = 3*(s*s*t)*p1.y
            let part4 = 3*(t*t*s)*p2.y
            let part44 = t*t*t*p3.y
            
            res.x = part1 + part11 + part2 + part22
            res.y = part3 + part33 + part4 + part44
            
            return res
        }
        
        func tangent(p0: PointF, p1: PointF, p2: PointF, p3: PointF, q0: PointF, q1: PointF) -> Double {
            var s: Double, r1: Double, r2: Double
            
            let A = cprod(p0: p0, p1: p1, p2: q0, p3: q1)
            let B = cprod(p0: p1, p1: p2, p2: q0, p3: q1)
            let C = cprod(p0: p2, p1: p3, p2: q0, p3: q1)
            
            let a = A - 2 * B + C
            let b = -2 * A + 2 * B
            let c = A
            
            let d = b * b - 4 * a * c
            
            if a == 0 || d < 0 {
                return -1.0
            }
            
            s = sqrt(d)
            
            r1 = (-b + s) / (2 * a)
            r2 = (-b - s) / (2 * a)
            
            if r1 >= 0 && r1 <= 1 {
                return r1
            } else if r2 >= 0 && r2 <= 1 {
                return r2
            } else {
                return -1.0
            }
        }
        
        func calcSums(path: Path) {
            var x: Int, y: Int
            path.x0 = path.pt[0].x
            path.y0 = path.pt[0].y
            
            path.sums = BoxedArray()
            let s = path.sums
            s.append(Sum(x: 0, y: 0, xy: 0, x2: 0, y2: 0))
            for i in 0 ..< path.len {
                x = path.pt[i].x - path.x0
                y = path.pt[i].y - path.y0
                s.append(Sum(x: s[i].x + x, y: s[i].y + y, xy: s[i].xy + x * y,
                             x2: s[i].x2 + x * x, y2: s[i].y2 + y * y))
            }
        }
        
        func calcLon(path: Path) {
            
            var n = path.len, pt = path.pt, dir: Int,
            pivk = [Int](repeating: 0, count: n),
            nc = [Int](repeating: 0, count: n),
            ct = [Int](repeating: 0, count: 4)
            
            path.lon = [Int](repeating: 0, count: n)
            
            var constraint = [PointI(x: 0, y: 0), PointI(x: 0, y: 0)],
            cur = PointI(x: 0, y: 0),
            off = PointI(x: 0, y: 0),
            dk = PointI(x: 0, y: 0),
            foundk: Int = 0
            
            var i: Int, j: Int, k1: Int, a: Int, b: Int, c: Int, d: Int, k: Int = 0
            
            for i in (0...n-1).reversed() {
                if pt[i].x != pt[k].x && pt[i].y != pt[k].y {
                    k = i + 1
                }
                nc[i] = k
            }
            
            for i in (0...n-1).reversed() {
                ct[0] = 0
                ct[1] = 0
                ct[2] = 0
                ct[3] = 0
                dir = (3 + 3 * (pt[mod(a: i + 1, n: n)].x - pt[i].x) +
                    (pt[mod(a: i + 1, n: n)].y - pt[i].y)) / 2
                ct[dir] += 1
                
                constraint[0].x = 0
                constraint[0].y = 0
                constraint[1].x = 0
                constraint[1].y = 0
                
                k = nc[i]
                k1 = i
                
                while true {
                    foundk = 0
                    dir =  (3 + 3 * sign(i: pt[k].x - pt[k1].x) +
                        sign(i: pt[k].y - pt[k1].y)) / 2
                    ct[dir] += 1
                    
                    if ct[0] != 0 && ct[1] != 0 && ct[2] != 0 && ct[3] != 0 {
                        pivk[i] = k1
                        foundk = 1
                        break
                    }
                    
                    cur.x = pt[k].x - pt[i].x
                    cur.y = pt[k].y - pt[i].y
                    
                    if xprod(p1: constraint[0], p2: cur) < 0 || xprod(p1: constraint[1], p2: cur) > 0 {
                        break
                    }
                    
                    if abs(cur.x) <= 1 && abs(cur.y) <= 1 {
                        // NOOP
                    } else {
                        off.x = cur.x + ((cur.y >= 0 && (cur.y > 0 || cur.x < 0)) ? 1 : -1)
                        off.y = cur.y + ((cur.x <= 0 && (cur.x < 0 || cur.y < 0)) ? 1 : -1)
                        if xprod(p1: constraint[0], p2: off) >= 0 {
                            constraint[0].x = off.x
                            constraint[0].y = off.y
                        }
                        off.x = cur.x + ((cur.y <= 0 && (cur.y < 0 || cur.x < 0)) ? 1 : -1)
                        off.y = cur.y + ((cur.x >= 0 && (cur.x > 0 || cur.y < 0)) ? 1 : -1)
                        if xprod(p1: constraint[1], p2: off) <= 0 {
                            constraint[1].x = off.x
                            constraint[1].y = off.y
                        }
                    }
                    k1 = k
                    k = nc[k1]
                    if (!cyclic(a: k, b: i, c: k1)) {
                        break
                    }
                }
                if foundk == 0 {
                    dk.x = sign(i: pt[k].x-pt[k1].x)
                    dk.y = sign(i: pt[k].y-pt[k1].y)
                    cur.x = pt[k1].x - pt[i].x
                    cur.y = pt[k1].y - pt[i].y
                    
                    a = xprod(p1: constraint[0], p2: cur)
                    b = xprod(p1: constraint[0], p2: dk)
                    c = xprod(p1: constraint[1], p2: cur)
                    d = xprod(p1: constraint[1], p2: dk)
                    
                    j = 10000000
                    if b < 0 {
                        j = Int(floor(Double(a / -b)))
                    }
                    if d > 0 {
                        j = min(j, Int(floor(Double(-c / d))))
                    }
                    pivk[i] = mod(a: k1+j, n: n)
                }
            }
            
            j=pivk[n-1]
            path.lon[n-1]=j
            for i in (0...n-2).reversed() {
                if cyclic(a: i+1, b: pivk[i], c: j) {
                    j=pivk[i]
                }
                path.lon[i]=j
            }
            
            i = n - 1
            while cyclic(a: mod(a: i+1, n: n), b: j, c: path.lon[i]) {
                i -= 1
            }
        }
        
        func bestPolygon(path: Path) {
            
            func penalty3(path: Path, i: Int, jj: Int) -> Double {
                
                let n = path.len, pt = path.pt
                let sums = path.sums
                var x: Int, y: Int, xy: Int, x2: Int, y2: Int,
                k: Int, a: Double, b: Double, c: Double, s: Double,
                px: Double, py: Double, ex: Int, ey: Int,
                r: Int = 0
                
                var j = jj
                
                if j >= n {
                    j -= n
                    r = 1
                }
                
                if r == 0 {
                    x = sums[j+1].x - sums[i].x
                    y = sums[j+1].y - sums[i].y
                    x2 = sums[j+1].x2 - sums[i].x2
                    xy = sums[j+1].xy - sums[i].xy
                    y2 = sums[j+1].y2 - sums[i].y2
                    k = j+1 - i
                } else {
                    x = sums[j+1].x - sums[i].x + sums[n].x
                    y = sums[j+1].y - sums[i].y + sums[n].y
                    x2 = sums[j+1].x2 - sums[i].x2 + sums[n].x2
                    xy = sums[j+1].xy - sums[i].xy + sums[n].xy
                    y2 = sums[j+1].y2 - sums[i].y2 + sums[n].y2
                    k = j+1 - i + n
                }
                
                px = (Double(pt[i].x + pt[j].x) / 2.0) - Double(pt[0].x)
                py = (Double(pt[i].y + pt[j].y) / 2.0) - Double(pt[0].y)
                ey = (pt[j].x - pt[i].x)
                ex = -(pt[j].y - pt[i].y)
                
                let a1 = (Double(x2) - Double(2*x)*px)
                let a2 = Double(k) + px*px
                a = a1 / a2
                
                let b1 = (Double(xy) - Double(x)*py - Double(y)*px)
                let b2 = Double(k) + px*py
                b = b1 / b2
                
                let c1 = (Double(y2) - Double(2*y)*py)
                let c2 = Double(k) + py*py
                c = c1 / c2
                
                let s1 = Double(ex * ex) * a
                let s2 = Double(2 * ex * ey) * b
                let s3 = Double(ey * ey) * c
                s = s1 + s2 + s3
                
                return sqrt(s)
            }
            
            var i: Int, j: Int, m: Int,
            n = path.len,
            pen = [Double](repeating: 0.0, count: n + 1),
            prev = [Int](repeating: 0, count: n + 1),
            clip0 = [Int](repeating: 0, count: n),
            clip1 = [Int](repeating: 0, count: n + 1),
            seg0 = [Int](repeating: 0, count: n + 1),
            seg1 = [Int](repeating: 0, count: n + 1),
            thispen: Double, best: Double, c: Int
            
            for i in 0..<n {
                c = mod(a: path.lon[mod(a: i-1, n: n)]-1, n: n)
                if c == i {
                    c = mod(a: i+1, n: n)
                }
                if c < i {
                    clip0[i] = n
                } else {
                    clip0[i] = c
                }
            }
            
            j = 1
            for i in 0..<n {
                while j <= clip0[i] {
                    clip1[j] = i
                    j += 1
                }
            }
            
            i = 0
            j = 0
            while i < n {
                seg0[j] = i
                i = clip0[i]
                j += 1
            }
            seg0[j] = n
            m = j
            
            i = n
            j = m
            while j > 0 {
                seg1[j] = i
                i = clip1[i]
                j -= 1
            }
            seg1[0] = 0
            
            pen[0]=0
            j = 1
            while j <= m {
                for i in seg1[j]...seg0[j] {
                    best = -1
                    for k in (clip1[i]...seg0[j-1]).reversed() {
                        thispen = penalty3(path: path, i: k, jj: i) + pen[k]
                        if best < 0 || thispen < best {
                            prev[i] = k
                            best = thispen
                        }
                    }
                    pen[i] = best
                }
                j += 1
            }
            path.m = m
            path.po = [Int](repeating: 0, count: m)
            
            i = n
            j = m - 1
            while i > 0 {
                i = prev[i]
                path.po[j] = i
                j -= 1
            }
        }

        func adjustVertices(path: Path) {
            
            func pointslope(path: Path, ii: Int, jj: Int, ctr: PointF, dir: PointF) {
                
                let n = path.len
                let sums = path.sums
                var x: Int, y: Int, xy: Int, x2: Int, y2: Int,
                k: Int, a: Double, b: Double, c: Double, lambda2: Double,
                l: Double,
                r: Int = 0

                var i = ii
                var j = jj
                
                while j>=n {
                    j-=n
                    r+=1
                }
                while i>=n {
                    i-=n
                    r-=1
                }
                while j<0 {
                    j+=n
                    r-=1
                }
                while i<0 {
                    i+=n
                    r+=1
                }
                
                x = sums[j+1].x-sums[i].x+r*sums[n].x
                y = sums[j+1].y-sums[i].y+r*sums[n].y
                x2 = sums[j+1].x2-sums[i].x2+r*sums[n].x2
                xy = sums[j+1].xy-sums[i].xy+r*sums[n].xy
                y2 = sums[j+1].y2-sums[i].y2+r*sums[n].y2
                k = j+1-i+r*n
                
                ctr.x = Double(x/k)
                ctr.y = Double(y/k)
                
                a = Double((x2-x*x/k)/k)
                b = Double((xy-x*y/k)/k)
                c = Double((y2-y*y/k)/k)
                
                lambda2 = (a+c+sqrt((a-c)*(a-c)+4*b*b))/2
                
                a -= lambda2
                c -= lambda2
                
                if abs(a) >= abs(c) {
                    l = sqrt(a*a+b*b)
                    if l != 0 {
                        dir.x = -b/l
                        dir.y = a/l
                    }
                } else {
                    l = sqrt(c*c+b*b)
                    if l != 0 {
                        dir.x = -c/l
                        dir.y = b/l
                    }
                }
                if l == 0 {
                    dir.x = 0.0
                    dir.y = 0.0
                }
            }
            
            var m = path.m, po = path.po, n = path.len, pt = path.pt,
                x0 = path.x0, y0 = path.y0,
                ctr = [PointF](repeating: PointF(x: 0, y: 0), count: m),
                dir = [PointF](repeating: PointF(x: 0, y: 0), count: m),
                q = [Quad](repeating: Quad(), count: m),
                v = [Double](repeating: 0.0, count: m), d: Double, j: Int,
                s = PointI(x: 0, y: 0)
            
            path.curve = Curve(n: m)
            
            for i in 0..<m {
                j = po[mod(a: i+1, n: m)]
                j = mod(a: j-po[i], n: n)+po[i]
                ctr[i] = PointF(x: 0, y: 0)
                dir[i] = PointF(x: 0, y: 0)
                pointslope(path: path, ii: po[i], jj: j, ctr: ctr[i], dir: dir[i])
            }
            
            for i in 0..<m {
                q[i] = Quad()
                d = dir[i].x * dir[i].x + dir[i].y * dir[i].y
                if d == 0.0 {
                    for j in 0..<3 {
                        for k in 0..<3 {
                            q[i].data[j * 3 + k] = 0
                        }
                    }
                } else {
                    v[0] = dir[i].y
                    v[1] = -dir[i].x
                    v[2] = -v[1] * ctr[i].y - v[0] * ctr[i].x
                    for l in 0..<3 {
                        for k in 0..<3 {
                            q[i].data[l * 3 + k] = v[l] * v[k] / d
                        }
                    }
                }
            }
            
            var Q: Quad, w: PointF, dx: Double, dy: Double,
                det: Double, min: Double, cand: Double, xmin: Double, ymin: Double
            
            for i in 0..<m {
                Q = Quad()
                w = PointF(x: 0, y: 0)
                
                s.x = pt[po[i]].x-x0
                s.y = pt[po[i]].y-y0
                
                j = mod(a: i-1, n: m)
                
                for l in 0..<3 {
                    for k in 0..<3 {
                        Q.data[l * 3 + k] = q[j].at(x: l, y: k) + q[i].at(x: l, y: k)
                    }
                }
                
                while true {
                    
                    det = Q.at(x: 0, y: 0)*Q.at(x: 1, y: 1) - Q.at(x: 0, y: 1)*Q.at(x: 1, y: 0)
                    if det != 0.0 {
                        w.x = (-Q.at(x: 0, y: 2)*Q.at(x: 1, y: 1) + Q.at(x: 1, y: 2)*Q.at(x: 0, y: 1)) / det
                        w.y = ( Q.at(x: 0, y: 2)*Q.at(x: 1, y: 0) - Q.at(x: 1, y: 2)*Q.at(x: 0, y: 0)) / det
                        break
                    }
                    
                    if (Q.at(x: 0, y: 0) > Q.at(x: 1, y: 1)) {
                        v[0] = -Q.at(x: 0, y: 1)
                        v[1] = Q.at(x: 0, y: 0)
                    } else if (Q.at(x: 1, y: 1) != 0) {
                        v[0] = -Q.at(x: 1, y: 1)
                        v[1] = Q.at(x: 1, y: 0)
                    } else {
                        v[0] = 1
                        v[1] = 0
                    }
                    d = v[0] * v[0] + v[1] * v[1]
                    v[2] = -v[1] * Double(s.y) - v[0] * Double(s.x)
                    for l in 0..<3 {
                        for k in 0..<3 {
                            Q.data[l * 3 + k] += v[l] * v[k] / d
                        }
                    }
                }
                dx = abs(w.x-Double(s.x))
                dy = abs(w.y-Double(s.y))
                if dx <= 0.5 && dy <= 0.5 {
                    path.curve.vertex[i] = PointF(x: w.x+Double(x0), y: w.y+Double(y0))
                    continue
                }
                
                min = quadform(Q: Q, w: PointI.toPointF(point: s))
                xmin = Double(s.x)
                ymin = Double(s.y)
                
                if (Q.at(x: 0, y: 0) != 0.0) {
                    for z in 0..<2 {
                        w.y = Double(s.y)-0.5+Double(z)
                        w.x = -(Q.at(x: 0, y: 1) * w.y + Q.at(x: 0, y: 2)) / Q.at(x: 0, y: 0)
                        dx = abs(w.x-Double(s.x))
                        cand = quadform(Q: Q, w: w)
                        if dx <= 0.5 && cand < min {
                            min = cand
                            xmin = w.x
                            ymin = w.y
                        }
                    }
                }
                
                if (Q.at(x: 1, y: 1) != 0.0) {
                    for z in 0..<2 {
                        w.x = Double(s.x)-0.5+Double(z)
                        w.y = -(Q.at(x: 1, y: 0) * w.x + Q.at(x: 1, y: 2)) / Q.at(x: 1, y: 1)
                        dy = abs(w.y-Double(s.y))
                        cand = quadform(Q: Q, w: w)
                        if dy <= 0.5 && cand < min {
                            min = cand
                            xmin = w.x
                            ymin = w.y
                        }
                    }
                }
                
                for l in 0..<3 {
                    for k in 0..<3 {
                        w.x = Double(s.x)-0.5+Double(l)
                        w.y = Double(s.y)-0.5+Double(k)
                        cand = quadform(Q: Q, w: w)
                        if cand < min {
                            min = cand
                            xmin = w.x
                            ymin = w.y
                        }
                    }
                }
                
                path.curve.vertex[i] = PointF(x: xmin + Double(x0), y: ymin + Double(y0))
            }
        }
        
        func reverse(path: Path) {
            var curve = path.curve, m = curve.n, v = curve.vertex,
            i: Int, j: Int, tmp: PointF
            
            i = 0
            j = m - 1
            while i < j {
                tmp = v[i]
                v[i] = v[j]
                v[j] = tmp
                i += 1
                j -= 1
            }
        }
        
        func smooth(path: Path) {
            let m = path.curve.n, curve = path.curve
            
            var j: Int, k: Int, dd: Double, denom: Double, alpha: Double,
            p2: PointF, p3: PointF, p4: PointF
            
            for i in 0..<m {

                j = mod(a: i+1, n: m)
                k = mod(a: i+2, n: m)
                p4 = interval(lambda: 1/2.0, a: curve.vertex[k], b: curve.vertex[j])
                
                denom = ddenom(p0: curve.vertex[i], p2: curve.vertex[k])
                if denom != 0.0 {
                    dd = dpara(p0: curve.vertex[i], p1: curve.vertex[j], p2: curve.vertex[k]) / denom
                    dd = abs(dd)
                    alpha = dd>1 ? (1 - 1.0/dd) : 0
                    alpha /= 0.75
                } else {
                    alpha = 4/3.0
                }
                curve.alpha0[j] = alpha
                
                if alpha >= info.alphamax {
                    curve.tag[j] = "CORNER"
                    curve.c[3 * j + 1] = curve.vertex[j]
                    curve.c[3 * j + 2] = p4
                } else {
                    if alpha < 0.55 {
                        alpha = 0.55
                    } else if alpha > 1 {
                        alpha = 1
                    }
                    p2 = interval(lambda: 0.5+0.5*alpha, a: curve.vertex[i], b: curve.vertex[j])
                    p3 = interval(lambda: 0.5+0.5*alpha, a: curve.vertex[k], b: curve.vertex[j])
                    curve.tag[j] = "CURVE"
                    curve.c[3 * j + 0] = p2
                    curve.c[3 * j + 1] = p3
                    curve.c[3 * j + 2] = p4
                }
                curve.alpha[j] = alpha
                curve.beta[j] = 0.5
            }
            curve.alphacurve = 1
        }
        
        func optiCurve(path: Path) {
            
            class Opti {
                var pen: Double
                var c: BoxedArray<PointF>
                var t: Double
                var s: Double
                var alpha: Double
                
                init() {
                    self.pen = 0
                    self.c = BoxedArray(repeating: PointF(x: 0, y: 0), count: 2)
                    self.t = 0
                    self.s = 0
                    self.alpha = 0
                }
            }
            
            func opti_penalty(path: Path, i: Int, j: Int, res: Opti,
                              opttolerance: Double, convc: Array<Int>, areac: Array<Double>) -> Int {
                
                var m = path.curve.n, curve = path.curve, vertex = curve.vertex,
                k: Int, k1: Int, k2: Int, conv: Int, i1: Int,
                area: Double, alpha: Double, d: Double, d1: Double, d2: Double,
                p0: PointF, p1: PointF, p2: PointF, p3: PointF, pt: PointF,
                A: Double, R: Double, A1: Double, A2: Double, A3: Double, A4: Double,
                s: Double, t: Double
                
                if i == j {
                    return 1
                }
                
                k = i
                i1 = mod(a: i+1, n: m)
                k1 = mod(a: k+1, n: m)
                conv = convc[k1]
                if (conv == 0) {
                    return 1
                }
                d = ddist(p: vertex[i], q: vertex[i1])

                k = k1
                while k != j {
                    k1 = mod(a: k+1, n: m)
                    k2 = mod(a: k+2, n: m)
                    if (convc[k1] != conv) {
                        return 1
                    }
                    if (sign(i: Int(cprod(p0: vertex[i], p1: vertex[i1], p2: vertex[k1], p3: vertex[k2]))) !=
                        conv) {
                        return 1
                    }
                    if (iprod1(p0: vertex[i], p1: vertex[i1], p2: vertex[k1], p3: vertex[k2]) <
                        d * ddist(p: vertex[k1], q: vertex[k2]) * -0.999847695156) {
                        return 1
                    }
                    
                    k = k1
                }
                
                p0 = curve.c[mod(a: i,n: m) * 3 + 2].copy()
                p1 = vertex[mod(a: i+1,n: m)].copy()
                p2 = vertex[mod(a: j,n: m)].copy()
                p3 = curve.c[mod(a: j,n: m) * 3 + 2].copy()
                
                area = areac[j] - areac[i]
                area -= dpara(p0: vertex[0], p1: curve.c[i * 3 + 2], p2: curve.c[j * 3 + 2])/2
                if (i>=j) {
                    area += areac[m]
                }
                
                A1 = dpara(p0: p0, p1: p1, p2: p2)
                A2 = dpara(p0: p0, p1: p1, p2: p3)
                A3 = dpara(p0: p0, p1: p2, p2: p3)
                
                A4 = A1+A3-A2    
                
                if (A2 == A1) {
                    return 1
                }
                
                t = A3/(A3-A4)
                s = A2/(A2-A1)
                A = A2 * t / 2.0
                
                if (A == 0.0) {
                    return 1
                }
                
                R = area / A
                alpha = 2 - sqrt(4 - R / 0.3)
                
                res.c[0] = interval(lambda: t * alpha, a: p0, b: p1)
                res.c[1] = interval(lambda: s * alpha, a: p3, b: p2)
                res.alpha = alpha
                res.t = t
                res.s = s
                
                p1 = res.c[0].copy()
                p2 = res.c[1].copy() 
                
                res.pen = 0
                
                k = mod(a: i+1,n: m)
                while k != j {
                    k1 = mod(a: k+1,n: m)
                    t = tangent(p0: p0, p1: p1, p2: p2, p3: p3, q0: vertex[k], q1: vertex[k1])
                    if (t < -0.5) {
                        return 1
                    }
                    pt = bezier(t: t, p0: p0, p1: p1, p2: p2, p3: p3)
                    d = ddist(p: vertex[k], q: vertex[k1])
                    if (d == 0.0) {
                        return 1
                    }
                    d1 = dpara(p0: vertex[k], p1: vertex[k1], p2: pt) / d
                    if (abs(d1) > opttolerance) {
                        return 1
                    }
                    if (iprod(p0: vertex[k], p1: vertex[k1], p2: pt) < 0 ||
                        iprod(p0: vertex[k1], p1: vertex[k], p2: pt) < 0) {
                        return 1
                    }
                    res.pen += d1 * d1
                    k=k1
                }
                
                k = i
                while k == j {
                    k1 = mod(a: k+1,n: m)
                    t = tangent(p0: p0, p1: p1, p2: p2, p3: p3, q0: curve.c[k * 3 + 2], q1: curve.c[k1 * 3 + 2])
                    if (t < -0.5) {
                        return 1
                    }
                    pt = bezier(t: t, p0: p0, p1: p1, p2: p2, p3: p3)
                    d = ddist(p: curve.c[k * 3 + 2], q: curve.c[k1 * 3 + 2])
                    if (d == 0.0) {
                        return 1
                    }
                    d1 = dpara(p0: curve.c[k * 3 + 2], p1: curve.c[k1 * 3 + 2], p2: pt) / d
                    d2 = dpara(p0: curve.c[k * 3 + 2], p1: curve.c[k1 * 3 + 2], p2: vertex[k1]) / d
                    d2 *= 0.75 * curve.alpha[k1]
                    if (d2 < 0) {
                        d1 = -d1
                        d2 = -d2
                    }
                    if (d1 < d2 - opttolerance) {
                        return 1
                    }
                    if (d1 < d2) {
                        res.pen += (d1 - d2) * (d1 - d2)
                    }
                    k=k1
                }
                
                return 0
            }
            
            var curve = path.curve, m = curve.n, vert = curve.vertex, 
            pt = Array<Int>(repeating: 0, count: m + 1),
            pen = Array<Double>(repeating: 0, count: m + 1),
            len = Array<Int>(repeating: 0, count: m + 1),
            opt = Array<Opti>(repeating: Opti(), count: m + 1),
            om: Int, i: Int,j: Int,r: Int,
            o = Opti(), p0: PointF,
            i1: Int, area: Double, alpha: Double, ocurve: Curve,
            s: Array<Double>, t: Array<Double>
            
            var convc = Array<Int>(repeating: 0, count: m),
                areac = Array<Double>(repeating: 0, count: m+1)
            
            for i in 0..<m {
                if (curve.tag[i] == "CURVE") {
                    let val = dpara(p0: vert[mod(a: i-1,n: m)], p1: vert[i], p2: vert[mod(a: i+1,n: m)])
                    convc[i] = sign(i: Int(val))
                } else {
                    convc[i] = 0
                }
            }
            
            area = 0.0
            areac[0] = 0.0
            p0 = curve.vertex[0]
            for i in 0..<m {
                i1 = mod(a: i+1, n: m)
                if (curve.tag[i1] == "CURVE") {
                    alpha = curve.alpha[i1]
                    area += 0.3 * alpha * (4-alpha) *
                        dpara(p0: curve.c[i * 3 + 2], p1: vert[i1], p2: curve.c[i1 * 3 + 2])/2
                    area += dpara(p0: p0, p1: curve.c[i * 3 + 2], p2: curve.c[i1 * 3 + 2])/2
                }
                areac[i+1] = area
            }
            
            pt[0] = -1
            pen[0] = 0
            len[0] = 0
            
            
            for j in 1...m {
                pt[j] = j-1
                pen[j] = pen[j-1]
                len[j] = len[j-1]+1
                
                i = j - 2
                while i>=0 {
                    r = opti_penalty(path: path, i: i, j: mod(a: j,n: m), res: o, opttolerance: info.opttolerance, convc: convc,
                                     areac: areac)
                    if r == 1 {
                        break
                    }
                    if (len[j] > len[i]+1 ||
                        (len[j] == len[i]+1 && pen[j] > pen[i] + o.pen)) {
                        pt[j] = i
                        pen[j] = pen[i] + o.pen
                        len[j] = len[i] + 1
                        opt[j] = o
                        o = Opti()
                    }
                    i -= 1
                }
            }
            om = len[m]
            ocurve = Curve(n: om)
            s = Array<Double>(repeating: 0, count: om)
            t = Array<Double>(repeating: 0, count: om)
            
            j = m
            for i in (0...om-1).reversed() {
                if (pt[j]==j-1) {
                    ocurve.tag[i]     = curve.tag[mod(a: j,n: m)]
                    ocurve.c[i * 3 + 0]    = curve.c[mod(a: j,n: m) * 3 + 0]
                    ocurve.c[i * 3 + 1]    = curve.c[mod(a: j,n: m) * 3 + 1]
                    ocurve.c[i * 3 + 2]    = curve.c[mod(a: j,n: m) * 3 + 2]
                    ocurve.vertex[i]  = curve.vertex[mod(a: j,n: m)]
                    ocurve.alpha[i]   = curve.alpha[mod(a: j,n: m)]
                    ocurve.alpha0[i]  = curve.alpha0[mod(a: j,n: m)]
                    ocurve.beta[i]    = curve.beta[mod(a: j,n: m)]
                    s[i] = 1.0
                    t[i] = 1.0
                } else {
                    ocurve.tag[i] = "CURVE"
                    ocurve.c[i * 3 + 0] = opt[j].c[0]
                    ocurve.c[i * 3 + 1] = opt[j].c[1]
                    ocurve.c[i * 3 + 2] = curve.c[mod(a: j,n: m) * 3 + 2]
                    ocurve.vertex[i] = interval(lambda: opt[j].s, a: curve.c[mod(a: j,n: m) * 3 + 2],
                                                b: vert[mod(a: j,n: m)])
                    ocurve.alpha[i] = opt[j].alpha
                    ocurve.alpha0[i] = opt[j].alpha
                    s[i] = opt[j].s
                    t[i] = opt[j].t
                }
                j = pt[j]
            }
            
            for i in 0..<om {
                i1 = mod(a: i+1,n: om)
                ocurve.beta[i] = s[i] / (s[i] + t[i1])
            }
            ocurve.alphacurve = 1
            path.curve = ocurve
        }
        
        for i in 0..<pathlist.count {
            let path = pathlist[i]
            calcSums(path: path)
            calcLon(path: path)
            bestPolygon(path: path)
            adjustVertices(path: path)
            
            if (path.sign == "-") {
                reverse(path: path)
            }
            
            smooth(path: path)
            
            if (info.optcurve) {
                optiCurve(path: path)
            }
        }
    }

    func getBezierPath(scale size: Double = 1.0) -> UIBezierPath {
        
        func path(curve: Curve, bezierPath: UIBezierPath) {
        
            func bezier(i: Int) {
                bezierPath.addCurve(to: CGPoint(x: (curve.c[i * 3 + 2].x * size), y: (curve.c[i * 3 + 2].y * size)),
                                    controlPoint1: CGPoint(x: (curve.c[i * 3 + 0].x * size), y: (curve.c[i * 3 + 0].y * size)),
                                    controlPoint2: CGPoint(x: (curve.c[i * 3 + 1].x * size), y: (curve.c[i * 3 + 1].y * size)))
            }

            func segment(i: Int) {
                bezierPath.addLine(to: CGPoint(x: (curve.c[i * 3 + 1].x * size), y: (curve.c[i * 3 + 1].y * size)))
                bezierPath.addLine(to: CGPoint(x: (curve.c[i * 3 + 2].x * size), y: (curve.c[i * 3 + 2].y * size)))
            }
            
            let n = curve.n
            
            let x = (curve.c[(n - 1) * 3 + 2].x * size)
            let y = (curve.c[(n - 1) * 3 + 2].y * size)
            bezierPath.move(to: CGPoint(x: x, y: y))

            for i in 0..<n {
                if (curve.tag[i] == "CURVE") {
                    bezier(i: i)
                } else if (curve.tag[i] == "CORNER") {
                    segment(i: i)
                }
            }
        }
        
        let len = pathlist.count
        var c: Curve
        
        let bezierPath = UIBezierPath()

        for i in 0 ..< len {
            c = pathlist[i].curve
            path(curve: c, bezierPath: bezierPath)
        }

        bezierPath.usesEvenOddFillRule = true
        return bezierPath
    }
    
    func getSVG(scale size: Double = 1.0, opt_type: String? = nil) -> String {
        
        func path(curve: Curve) -> String {
            
            func bezier(i: Int) -> String {
                var b = "C " + String(format: "%f", (curve.c[i * 3 + 0].x * size).toFixed(3)) + " " +
                    String(format: "%f",(curve.c[i * 3 + 0].y * size).toFixed(3)) + ","
                b += String(format: "%f", (curve.c[i * 3 + 1].x * size).toFixed(3)) + " " +
                    String(format: "%f", (curve.c[i * 3 + 1].y * size).toFixed(3)) + ","
                b += String(format: "%f", (curve.c[i * 3 + 2].x * size).toFixed(3)) + " " +
                    String(format: "%f", (curve.c[i * 3 + 2].y * size).toFixed(3)) + " "
                return b
            }
            
            func segment(i: Int) -> String {
                var s = "L " + String(format: "%f", (curve.c[i * 3 + 1].x * size).toFixed(3)) + " " +
                    String(format: "%f", (curve.c[i * 3 + 1].y * size).toFixed(3)) + " "
                s += String(format: "%f", (curve.c[i * 3 + 2].x * size).toFixed(3)) + " " +
                    String(format: "%f", (curve.c[i * 3 + 2].y * size).toFixed(3)) + " "
                return s
            }
            
            let n = curve.n
            var p = "M" + String(format: "%f", (curve.c[(n - 1) * 3 + 2].x * size).toFixed(3)) +
                " " + String(format: "%f", (curve.c[(n - 1) * 3 + 2].y * size).toFixed(3)) + " "
            for i in 0..<n {
                if (curve.tag[i] == "CURVE") {
                    p += bezier(i: i)
                } else if (curve.tag[i] == "CORNER") {
                    p += segment(i: i)
                }
            }

            return p
        }
        
        var w = Double(bm.w) * size, h = Double(bm.h) * size,
        len = pathlist.count, c: Curve, strokec: String, fillc: String, fillrule: String
        
        var svg = "<svg id=\"svg\" version=\"1.1\" width=\"\(w)\" height=\"\(h)\" xmlns=\"http://www.w3.org/2000/svg\">"
        svg += "<path d=\""
        for i in 0 ..< len {
            c = pathlist[i].curve
            svg += path(curve: c)
        }
        if (opt_type == "curve") {
            strokec = "black"
            fillc = "none"
            fillrule = ""
        } else {
            strokec = "none"
            fillc = "black"
            fillrule = " fill-rule=\"evenodd\""
        }
        svg += "\" stroke=\"\(strokec)\" fill=\"\(fillc)\" \(fillrule)/></svg>"
        return svg
    }
}

extension Double {
    func toFixed(_ places: Int) -> Double {
        let divisor = pow(10.0, Double(places))
        return (self * divisor).rounded() / divisor
    }
}

class BoxedArray<T> : MutableCollection {
    var array: [T]
    
    init() {
        self.array = [T]()
    }
    
    init(repeating: T, count: Int) {
        self.array = [T](repeating: repeating, count: count)
    }
    
    func index(after i: (Int)) -> (Int) {
        return array.index(after: i)
    }
    
    var endIndex: (Int) {
        return array.endIndex
    }
    
    var startIndex: (Int) {
        return array.startIndex
    }
    
    func append(_ item: T) {
        array.append(item)
    }
    
    subscript (index: Int) -> T {
        get { return array[index] }
        set(newValue) { array[index] = newValue }
    }
}


