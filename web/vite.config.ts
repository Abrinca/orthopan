import {defineConfig} from 'vite';

export default defineConfig({
  base: './',
  server: {
    proxy: {
      '/api': 'http://127.0.0.1:8000',
    },
  },
  build: {
    // Ship the UI inside the Python package so wheels serve it without bun.
    outDir: '../src/orthopan/static',
    emptyOutDir: true,
    sourcemap: true,
    chunkSizeWarningLimit: 5000,
  },
});
