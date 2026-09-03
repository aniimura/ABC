# [Falt1] トラックのゴール(2026-09-04 設定、ABC3c)

対象: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299
(物理46ページ、`papers.json`に短縮タグ`Falt1`で登記済み、`0_Source`にPDF/txtあり)。

**★構造は実測済み(pdftotext + 260dpi目視1点)だが、逐語は大部分まだ目視していない。**
JSTORスキャンのpdftotextは**地の文の単語間スペースが失われる**(例:
"We intend to describe" → "Weintendtodescribea")——数式・定理番号は無事だが、
散文コメントの逐語引用にはPDF目視が必須。260dpiで1ページ(物理p.17、
LocProPが直接引用するTheorem 4.4を含む)を確認し、定理本体は完全に読めることを確認した。

## 0. ゴール(現在地)

> **Falt1 Chapter I §1 0/2, §2 0/4, §3 0/2, §4 0/5,
> Chapter II §1 0/2, §2 0/4, §3 0/2,
> Chapter III §1 0/3, §2 0/2, §3 0/3, §4 0/1+ —— 合計 0/30+**

★Chapter III §4は未確認(ファイル末尾に近く、もう1-2件ある可能性)。
★Chapter番号(I/II/III)はpdftotextの見出し検出に失敗した(空白喪失)ため、
節番号が1に戻る境界(§4→§1)を手がかりに3つに分割した——**Chapter番号の
目視確認はまだ**。

## 1. 番号付け規則

★★**この論文は "N.M. Kind" の順**(番号が先、種別が後)。LocProP等の
"Kind N.M"(種別が先)と**逆**。`tools/paper-items.mjs`はそのまま使えない
(専用の抽出スクリプトを`scratchpad/falt1-items.mjs`に書いた、再現可能)。

## 2. LocProPが直接引用する箇所(確認済み)

| LocProP側 | Falt1側 | 状態 |
|---|---|---|
| [LocProP] Lemma 2.1 | Chapter I, Theorem 4.4, (i)(iv) | ★★物理p.17(印字p.270)で目視確認済み。式は完全に読める |
| [LocProP] Lemma 2.3 | Chapter I(Leray-Serre + 前2つ) | 未確認 |
| [LocProP] 全体(§2の土台) | almost étale extension理論(§I-2「Almost unramified extensions」) | 未確認 |

## 3. 次にやること

1. 残り約45ページの260dpi目視(散文部分の逐語確認)。特に Chapter I §2
   (almost étale extensionの定義そのもの、LocProPの土台)が最優先。
2. Chapter番号(I/II/III)を見出しから直接確認する(現在は節番号の巡回から推定)。
3. Chapter III §4の末尾を確認し、正確な合計項目数を確定する。
4. Skeleton/Falt1/ChapterI.lean 等を作り、Theorem 4.4(LocProPが直接消費する)
   から着手する。

関連: [[measure-mathlib-before-skeleton]] / `ResearchPaper/mathlib-gap.json`の
`locprop-perfectoid-hodge-tate`(この論文が主典拠)
