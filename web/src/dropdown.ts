// Position-fixed dropdowns anchored to a control in the (horizontally
// scrolling) top bar: place below the anchor, clamped to the viewport.

export function placeDropdown(anchor: Element, el: HTMLElement): void {
  const r = anchor.getBoundingClientRect();
  // constrain size first (overrides CSS min-width on narrow viewports), then measure
  const maxW = Math.max(0, window.innerWidth - 8);
  el.style.maxWidth = `${maxW}px`;
  el.style.minWidth = `min(${maxW}px, var(--dropdown-min, 0px))`;
  el.style.maxHeight = `${Math.max(80, window.innerHeight - r.bottom - 8)}px`;
  const w = el.offsetWidth || 230;
  const left = Math.max(4, Math.min(Math.round(r.left), window.innerWidth - w - 4));
  el.style.left = `${left}px`;
  el.style.top = `${Math.round(r.bottom + 2)}px`;
}

/** Re-place `el` under `anchor` whenever `scroller` scrolls or the window resizes, while `isOpen()`. */
export function trackDropdown(anchor: Element, el: HTMLElement, scroller: Element, isOpen: () => boolean): void {
  const update = () => {
    if (isOpen()) placeDropdown(anchor, el);
  };
  scroller.addEventListener('scroll', update, {passive: true});
  window.addEventListener('resize', update);
}
