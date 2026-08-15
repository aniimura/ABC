// 各論文の参考文献欄を解いて、引用キー → 0_Source のファイル名 を作る。
//
// ★なぜ要るか(2026-08-15 実測)
//   IUT の4本は [IUTchI] のような**記号的な**キーを使うが、古い論文は
//   [Mzk3] のような**数字型**のキーを使う。しかも**キーは論文ごとに独立**である
//   (SemiAnbd の [Mzk4] と EtTh の [Mzk4] は別物でありうる)。
//   我々の依存グラフは記号的キーしか解決していなかったので、
//   古い論文からの辺を**ごっそり落としていた**:
//     AbsTopIII 318 件 / EtTh 254 件 / SemiAnbd 92 件 / FrdI 58 件 / IUTchIII 0 件
//   ★実害: anabelioid の定義元([Mzk4] = The Geometry of Anabelioids)が
//    グラフに1度も現れなかった。
//
// ★限界: 参考文献欄の題名と 0_Source のファイル名を**正規化して照合**している。
//   題名が一致しないものは解決できない(望月氏以外の文献はそもそも 0_Source に無い)。
//   pdftotext 由来なので題名が壊れている可能性もある。**目視していない。**

import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { join } from 'node:path';

const norm = (s) => s.toLowerCase()
  .replace(/[éèê]/g, 'e').replace(/[üú]/g, 'u').replace(/ö/g, 'o')
  .replace(/[^a-z0-9]+/g, ' ').replace(/\s+/g, ' ').trim();

/** 0_Source の本体 PDF/TXT のファイル名(拡張子なし)一覧 */
export function sourceTitles(SRC) {
  return readdirSync(SRC)
    .filter((f) => f.endsWith('.txt'))
    .map((f) => f.slice(0, -4))
    .filter((f) => !/\(comments\)|\(related|LaTeX version/.test(f));
}

/**
 * 1論文の参考文献欄から key → 0_Source のファイル名 を作る。
 * 参考文献の行は `[Key] Author, Title, ...` の形。題名は最初のカンマ以降・次のカンマまで。
 */
export function bibOf(SRC, file, titles) {
  const p = join(SRC, `${file}.txt`);
  if (!existsSync(p)) return new Map();
  const txt = readFileSync(p, 'utf8').split(/\r?\n/);
  const normTitles = titles.map((t) => [norm(t), t]);
  const out = new Map();
  for (let i = 0; i < txt.length; i++) {
    const m = txt[i].match(/^\s*\[([A-Za-z][A-Za-z0-9]*)\]\s+(.*)$/);
    if (!m) continue;
    // 題名は次の `[` 行まで(折り返しがあるため数行つなぐ)
    let body = m[2];
    for (let j = i + 1; j < Math.min(i + 5, txt.length); j++) {
      if (/^\s*\[[A-Za-z]/.test(txt[j])) break;
      body += ' ' + txt[j];
    }
    // ★題名の切り出し: 書式は `[Key] Author, Title, Publisher/Journal …`。
    //   最初のカンマまでが著者。以降の最初の「, 」までを題名とみなす。
    const after = body.replace(/^[^,]*,\s*/, '');
    const title = after.split(/,\s/)[0];
    const nt0 = norm(title);
    if (nt0.length < 10) continue;

    // ★部分文字列一致は弱すぎる(2026-08-15 実測)。実害:
    //   [Mzk3] "The Absolute Anabelian Geometry of Hyperbolic Curves"(0_Source に無い)が
    //     "Absolute Anabelian Geometry" に誤マッチ
    //   [Mzk8] "Galois Sections in Absolute Anabelian Geometry" も同じ題名に誤マッチ
    //   → **接頭辞であること**と**題名の 60% 以上を覆うこと**を要求する。
    //     満たさないものは解決しない。★誤った辺は欠けた辺より悪い。
    const strip = (s) => s.replace(/^the /, '');
    const nb = strip(nt0);
    let best = null;
    for (const [nt, t] of normTitles) {
      const c = strip(nt);
      if (c.length < 10) continue;
      if (!nb.startsWith(c)) continue;
      if (c.length / nb.length < 0.6) continue;
      if (!best || c.length > strip(norm(best)).length) best = t;
    }
    if (best && !out.has(m[1])) out.set(m[1], best);
    else if (!best) { if (!out.has(`?${m[1]}`)) out.set(`?${m[1]}`, `UNRESOLVED:${title.slice(0, 60)}`); }
  }
  return out;
}

/** 全論文ぶんの bib を作る。tagOfFile は ファイル名 → 我々のタグ。 */
export function buildBibs(SRC, fileOf) {
  const titles = sourceTitles(SRC);
  const tagOfFile = new Map(Object.entries(fileOf).map(([tag, f]) => [f, tag]));
  const bibs = new Map();   // tag -> Map(key -> tag)
  for (const [tag, file] of Object.entries(fileOf)) {
    const raw = bibOf(SRC, file, titles);
    const m = new Map();
    for (const [k, title] of raw) {
      const t = tagOfFile.get(title);
      m.set(k, t ?? `FILE:${title}`);   // 我々のタグが無い論文は FILE: 付きで返す
    }
    bibs.set(tag, m);
  }
  return bibs;
}
