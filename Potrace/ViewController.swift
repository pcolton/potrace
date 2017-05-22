//
//  ViewController.swift
//  Potrace
//
//  Created by Paul Colton on 5/19/17.
//  Copyright © 2017 Paul Colton. All rights reserved.
//

import UIKit

class ViewController: UIViewController {
    @IBOutlet weak var original: UIImageView!
    @IBOutlet weak var curves: UIImageView!

    override func viewDidLoad() {
        super.viewDidLoad()

        if let originalImage = original.image,
            let originalImageData = originalImage.pixelData() {
            
            let potrace = Potrace(data: UnsafeMutableRawPointer(mutating: originalImageData),
                                  width: Int(originalImage.size.width),
                                  height: Int(originalImage.size.height))
            
            potrace.process()
            
            let bezier = potrace.getBezierPath(scale: 2.0)
            let newImage = shapeImageWithCGPath(path: bezier.cgPath,
                                                size: curves.frame.size,
                                                fillColor: nil,
                                                strokeColor: UIColor.black)
            
            
            curves.image = newImage
        }
    }

    func shapeImageWithCGPath(path: CGPath, size: CGSize, fillColor: UIColor?, strokeColor: UIColor?, strokeWidth: CGFloat = 0.0) -> UIImage! {
        
        UIGraphicsBeginImageContext(size)
        let context = UIGraphicsGetCurrentContext()
        var image = UIImage()
        if let context  = context {
            context.saveGState()
            context.addPath(path)
            if strokeColor != nil {
                strokeColor!.setStroke()
                context.setLineWidth(strokeWidth)
            } else { UIColor.clear.setStroke() }
            fillColor?.setFill()
            context.drawPath(using: .fillStroke)
            image = UIGraphicsGetImageFromCurrentImageContext()!
            context.restoreGState()
            UIGraphicsEndImageContext()
        }
        return image
    }

}

extension UIImage {
    func pixelData() -> [UInt8]? {
        let size = self.size
        let dataSize = size.width * size.height * 4
        var pixelData = [UInt8](repeating: 0, count: Int(dataSize))
        let colorSpace = CGColorSpaceCreateDeviceRGB()
        let context = CGContext(data: &pixelData,
                                width: Int(size.width),
                                height: Int(size.height),
                                bitsPerComponent: 8,
                                bytesPerRow: 4 * Int(size.width),
                                space: colorSpace,
                                bitmapInfo: CGImageAlphaInfo.noneSkipLast.rawValue)
        guard let cgImage = self.cgImage else { return nil }
        context?.draw(cgImage, in: CGRect(x: 0, y: 0, width: size.width, height: size.height))
        
        return pixelData
    }
}


