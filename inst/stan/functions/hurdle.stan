vector hurdle_smooth(vector x, real threshold, real sharpness) {
  return((threshold - x) / threshold * sharpness);
}
