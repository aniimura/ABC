/**
 * tools/heading-shape.mjs —— 「原典の `.txt` で、その項目の見出しはどの行か」の**唯一の実装**
 * ==========================================================================
 *
 * ★なぜ別ファイルに出したか（2026-09-07 メタ第 20 回、backlog M80）
 * ------------------------------------------------------------------
 * この判定（M59 の `headingLineShape` / M65 の `(?![a-z])` / M71 の「形が合う最初の行を採る」）は
 * これまで `tools/brief.mjs` の**内部関数**だった。M76 が見張りを足したが、
 * ★見張りは `brief.mjs --audit-proof` の側にしか無く、**ゲートで落ちない**。
 * ★★落ちない見張りは、やがて誰も見なくなる（本体の判断。G7 の基準・G1 の `Found/` 非対称と同じ話）。
 *
 * `check.mjs` の selftest から叩けるようにするには 2 つしか道が無い:
 *   (a) `check.mjs` に判定を**写す** … ★写した方を測ることになる。**採らない**
 *       （`--audit-proof` を `brief.mjs` に相乗りさせた理由と同じ。M47-1 のコメント参照）。
 *   (b) 判定を**共有モジュールに 1 本化**し、両方が同じ実装を叩く … ★こちら。
 * ★`brief.mjs` を直接 import することはできない（モジュール本体が CLI で、
 *   引数が無いと `process.exit(2)` する）。★`entities.mjs`（M55）と同じ切り出し方である。
 *
 * ★★`check.mjs` の `SELF_HASH` には**入れない**（意図的）。
 *   鍵の不変条件は「`squash()`（pdftotext 出力の正規化）に効くものは全部鍵に入っている」であり、
 *   ★このファイルは `.txt` の見出し判定にしか効かず、PDF キャッシュの中身には 1 バイトも効かない。
 *   （`entities.mjs` は `decodeEntities` 経由で `squash` に効くので鍵に入っている。★対称ではない。）
 */

const HEADWORDS = 'Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example|Claim|Fact|Exercise';

/** 原典の見出し行を捕まえる正規表現（`.txt` は行頭に見出しが立つ——48 本で実測）。 */
export function headingRe() {
  // ★★M51(メタ第 14 回): Faltings は見出しを**番号先行**で書く(`1.2. Theorem.`)。
  //   旧式は `<種別> <番号>` しか見ておらず、Falt1 の 13 項目が**無音で失敗**していた。
  //   ★他の 15 論文は 348/360 が outcome も抽出行数も完全一致(副作用 0)。
  return new RegExp('^(?:(' + HEADWORDS + ')\\s+[0-9]+(?:\\.[0-9]+)*\\.?'
    + '|[0-9]+(?:\\.[0-9]+)*\\.\\s*(?:' + HEADWORDS + ')\\b)');
}

/**
 * ★**「見出しに見える行」が本当に見出しかを、行の形だけで見分ける**（メタ第 16 回 M59）。
 *
 * 動機（実測 2026-09-07）: `proofParagraphOf` が「終端記号に当たらず切った」54 件は
 * **全部が境界で切っており、うち 36 件は文の途中**だった。切った行の実物:
 *
 *     Proposition 4.4(ii), which shows v(β) = qi.      ← [Yoshida08] Prop 6.14 の段 1
 *     Definition 1.3, (v), (b), or the essential …     ← [FrdI] Prop 1.4
 *     Remark 1.1.1; the factorization of co-angular …  ← [FrdI] Prop 1.10（×6 項目）
 *
 * ★これらは**証明の中の相互参照が行頭に折り返しただけ**である。既存の
 * 「同じ鍵の最初の出現だけを境界にする」規則（M47-2）が効かなかったのは、
 * ★**鍵が末尾のドットで割れている**ため —— 真の見出しは `Proposition 4.4.`、
 * 相互参照は `Proposition 4.4`（直後が `(` や `,`）で、**別の鍵として両方が境界になっていた**
 * （[FrdI] は 115 鍵中 32 鍵がこの対で、対の片割れは全部が相互参照）。
 *
 * 見分けは**行の形だけ**でつく（前の行は見ない。走り込みの導入文
 * 「Then we make the following」があるので前の行は当てにならない）:
 *   - 番号先行（`2.2. Theorem.…` = [Falt1]）… 相互参照はこの形を取らない ⇒ 見出し
 *   - `Proposition 1.2:` … [pGC] の書式 ⇒ 見出し
 *   - `Definition 1.1.` のように**ドットで閉じ**、後ろが空か大文字か `(` ⇒ 見出し
 *   - `Proposition 6.6 (Sen [14]).` のように**空白 + `(`** ⇒ 見出し
 *   - それ以外（`, (v)` `; the` `(ii),` ` that` が続く）⇒ **相互参照。境界にしない**
 *
 * ★実測: 誤った境界 36/36 を落とし、**真の見出しを失う鍵は 1 個だけ**
 * （[Tate] `Theorem 1`。★この論文の `.txt` は 2 段組の OCR で、
 * そもそも真の見出しが 1 行も残っていない —— M39 測定 5 で既知）。
 */
export function headingLineShape(line, m) {
  const hit = m[0];
  if (/^[0-9]/.test(hit)) return true;             // 番号先行（[Falt1]）
  const rest = line.slice(hit.length);
  if (/^\s*:/.test(rest)) return true;             // [pGC] `Proposition 1.2:`
  if (hit.endsWith('.')) return rest.trim() === '' || /^\s*[(“"A-Z]/.test(rest);
  return /^\s+\(/.test(rest);                      // `Proposition 6.6 (Sen [14]).`
}

/**
 * 項目（`kind num`）の見出しを探す 2 本の正規表現。
 *
 * ★★M62/M64(メタ第 17 回): 見出し語の直後は `\b` ではなく `(?![a-z])` で閉じる。
 *   [Falt1] の `.txt` 279 行は `2.1. DefinitionS.uppose A is a ring` ——
 *   ★**OCR がピリオドを 1 文字ずらしている**（`Definition.` `Suppose` → `DefinitionS` `.uppose`）ので、
 *   `Definition` の直後が語構成文字になり `\b` が立たず、**無音で `noheading` になっていた**。
 *   `(?![a-z])` なら `DefinitionS` を通し、`Definitions`（複数形＝見出しではない）は弾く。
 *
 * ★★落とし穴（M66 が予告していた）: `want` は `(?!\.?[0-9])` が**先読み**なので
 *   **末尾のドットを消費しない**。そのまま `headingLineShape` に掛けると `m[0]` が
 *   `Theorem 3.11`（ドット無し）になり、★**真の見出し `Theorem 3.11.` まで落ちる。**
 *   ⇒ 形の判定用に**ドットを取り込む** `wantDot` を別に作る。
 *   ★`\.?` を先読みより前に置いてはいけない（`Theorem 3.11.2` が backtrack で通る）。
 */
export function itemHeadingRes(kind, num) {
  const numEsc = String(num).replace(/\./g, '\\.');
  return {
    want: new RegExp(`^(?:${kind}\\s+${numEsc}(?!\\.?[0-9])|${numEsc}\\.\\s*${kind}(?![a-z]))`),
    wantDot: new RegExp(`^(?:${kind}\\s+${numEsc}(?!\\.?[0-9])\\.?|${numEsc}\\.\\s*${kind}(?![a-z]))`),
  };
}

/**
 * ★★★**「行頭に立ち、かつ見出しの形をしている行」を全部数え、最初のものを採る**（M71 + M76）。
 *
 * M71 の申告: 「最初に形が合う行」で当たるのは、★**原典が項目を昇順に一度だけ立てるから**にすぎない。
 *   [IUTchIII] 11870 行 `Theorem 3.11. For n ∈Z, write` は
 *   ★**折り返した相互参照なのに見出しの形をしている**（直前の 11869 行が
 *   `… we are in the situation of` で終わっている）。今回落ちなかったのは
 *   ★真の見出し（10489 行）が**先に**出るからにすぎない。
 *   ⇒ ★序文で先に形の合う引用が出る論文があれば破れる。
 *
 * ⇒ ★**形が合う行が 2 箇所以上ある項目を数える**（`hits`）。これは推測を含まない**厳密な数**であり、
 *   「破れうる持ち場」の**上界**である（破れるとしたら必ずこの集合の中で起きる）。
 *   `--audit-proof` が毎回この数と明細を出し、★`check.mjs --selftest` の `txtCases` が
 *   ★**規則そのもの（形の判定 ＋ 最初を採る）をゲートで毎回叩く**（メタ第 20 回 M80）。
 *   ★ヒューリスティックで「怪しい」と言わないのが要点
 *   （直前行の形で当てにいく案も測ったが、[CorrHyp]/[LocProP]/[pGC] の
 *    `Then we make the following` ＋ `Definition 1.1.` という**正しい**組版を
 *    10 件中 5 件で誤って挙げた ＝ 偽陽性 50%。★採らなかった）。
 *
 * @param {string[]} lines `.txt` を `\n` で割ったもの
 * @returns {{at:number, hits:number[]}} `at` は採用行（**1 始まり**。無ければ 0）、
 *   `hits` は形が合った行の**1 始まり**の一覧（昇順）。
 */
export function pickHeading(lines, kind, num) {
  const { wantDot } = itemHeadingRes(kind, num);
  const hits = [];
  for (let i = 0; i < lines.length; i++) {
    const m = wantDot.exec(lines[i]);
    if (m && headingLineShape(lines[i], m)) hits.push(i + 1);
  }
  // ★ここが M76 の「破れる形」そのもの。`hits[0]`（＝**最初**）を採る。
  //   `hits[hits.length - 1]` に変えると [IUTchIII] が 11/13 → 0/13 に落ちる（M76 が実測）。
  return { at: hits.length ? hits[0] : 0, hits };
}
