const puppeteer = require('puppeteer');
const path = require('path');
(async () => {
  const browser = await puppeteer.launch({
    headless: 'new',
    args: ['--no-sandbox', '--disable-setuid-sandbox', '--font-render-hinting=none']
  });
  const page = await browser.newPage();
  await page.setViewport({ width: 1600, height: 2000, deviceScaleFactor: 3 });
  const htmlPath = path.resolve(__dirname, 'figure1.html');
  await page.goto(`file://${htmlPath}`, { waitUntil: 'networkidle0' });
  await new Promise(r => setTimeout(r, 500));
  const bbox = await page.evaluate(() => {
    const el = document.getElementById('figure');
    const r = el.getBoundingClientRect();
    return { x: r.x, y: r.y, width: r.width, height: r.height };
  });
  const outDir = process.argv[2] || __dirname;
  await page.screenshot({
    path: path.join(outDir, 'Figure1_LSEA_algorithm.png'),
    clip: { x: bbox.x, y: bbox.y, width: bbox.width, height: bbox.height },
    omitBackground: false
  });
  await page.pdf({
    path: path.join(outDir, 'Figure1_LSEA_algorithm.pdf'),
    width: `${bbox.width}px`,
    height: `${bbox.height + 20}px`,
    printBackground: true,
    margin: { top: 0, right: 0, bottom: 0, left: 0 }
  });
  await browser.close();
  console.log(`Rendered: ${bbox.width}x${bbox.height}px (3x = ${bbox.width*3}x${bbox.height*3}px)`);
  console.log(`PNG: ${path.join(outDir, 'Figure1_LSEA_algorithm.png')}`);
  console.log(`PDF: ${path.join(outDir, 'Figure1_LSEA_algorithm.pdf')}`);
})();
