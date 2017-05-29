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
    @IBOutlet weak var curveTolerance: UILabel!
    @IBOutlet weak var turdSize: UILabel!
    @IBOutlet weak var turdSizeSlider: UISlider!
    @IBOutlet weak var curveToleranceSlider: UISlider!
    @IBOutlet weak var optCurveSwitch: UISwitch!
    
    var potrace: Potrace!
    
    var width: Int!
    var height: Int!
    var imagePixels: [UInt8]!

    override func viewDidLoad() {
        super.viewDidLoad()

        curveTolerance.text = String("Curve Tolerance (\(curveToleranceSlider.value))")
        turdSize.text = String("Turd Size (\(Int(round(turdSizeSlider.value))))")
        curveToleranceSlider.isEnabled = optCurveSwitch.isOn

        if let originalImage = original.image,
            let pixels = originalImage.pixelData() {

            self.width = Int(originalImage.size.width)
            self.height = Int(originalImage.size.height)
            self.imagePixels = pixels
            
            updateImage(settings: Potrace.Settings())
        }
    }
    
    func updateImage(settings: Potrace.Settings) {
        self.potrace = Potrace(data: UnsafeMutableRawPointer(mutating: self.imagePixels),
                               width: self.width,
                               height: self.height)
        
        self.potrace.process(settings: settings)

        let bezier = potrace.getBezierPath(scale: 2.0)
        
        DispatchQueue.main.async {
            let newImage = self.imageFromBezierPath(path: bezier, size: self.curves.frame.size)
            self.curves.image = newImage
        }
    }
    
    @IBAction func turdSizeTouchesUp(_ sender: UISlider) {
        updateImage(settings: collectSettings())
    }
    
    @IBAction func optCurveChanged(_ sender: UISwitch) {
        curveToleranceSlider.isEnabled = sender.isOn
        updateImage(settings: collectSettings())
    }

    @IBAction func turdSizeValueChanged(_ sender: UISlider) {
        turdSize.text = String("Turd Size (\(Int(sender.value))")
    }
    
    @IBAction func optToleranceSlider(_ sender: UISlider) {
        curveTolerance.text = String("Curve Tolerance (\(sender.value))")
    }

    @IBAction func optToleranceTouchesUp(_ sender: UISlider) {
        updateImage(settings: collectSettings())
    }

    func collectSettings() -> Potrace.Settings {
    
        var settings = Potrace.Settings()
        settings.turdsize = Int(turdSizeSlider.value)
        settings.optcurve = optCurveSwitch.isOn
        settings.opttolerance = Double(curveToleranceSlider.value)
        
        return settings
    }
    
    func imageFromBezierPath(path: UIBezierPath, size: CGSize) -> UIImage {
        var image = UIImage()
        UIGraphicsBeginImageContext(size)
        if let context = UIGraphicsGetCurrentContext() {
            context.saveGState()
            path.fill()
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
