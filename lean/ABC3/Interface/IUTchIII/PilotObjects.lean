import ABC3.Meta.Claim
import Mathlib.Data.Real.Basic
import Mathlib.Data.Set.Lattice
import Mathlib.Order.WithBot
import Mathlib.Topology.Compactness.Compact

/-!
# [IUTchIII] Corollary 3.12 — まだ構築できていない基礎

Corollary 3.12(物理 p.173–174)が**所与として読み出す**もの。
構造化: `ResearchPaper/1_Structured/Inter-universal Teichmuller Theory III/section-3.html`
(PDF 目視確認済み: 物理 p.153・154・156・173・174)。

`axiom` ではなく `structure` で受ける。

## ★このファイルの設計上いちばん大事なこと — 非対称性

原文は Θ 側と q 側を**違う扱い**にしている。逐語:

原文 (IUTchIII p.174):
> relative to the relevant Kummer isomorphisms [cf. Theorem 3.11, (ii)], in the
> multiradial representation of Theorem 3.11, (i), which we regard as subject to
> the indeterminacies (Ind1), (Ind2), (Ind3) described in Theorem 3.11, (i), (ii).

原文 (IUTchIII p.174):
> which we do not regard as subject to the indeterminacies (Ind1), (Ind2),
> (Ind3) described in Theorem 3.11, (i), (ii).

前者が Θ 側、後者が q 側である。**この非対称性が消えると不等式は内容を失う。**

非対称性は装飾ではなく、**取る対象そのものが違う**という形で原文に現れている:

- Θ 側 — 「the **union** of the **possible images** of a Θ-pilot object」の
  「holomorphic **hull**」の対数体積。像が**複数ありうる**。
- q 側 — 「the **image** of a q-pilot object」の対数体積。**単数**。

ゆえに型でこう写す:

```
possibleThetaImages : Set (Set Amb)   -- 原文「the possible images」= 像の集まり
qImage              : Set Amb         -- 原文「the image」= 単数
```

そして Θ 側は `logVol (hull (⋃₀ possibleThetaImages))`、q 側は `logVol qImage`。
**2つの実数 a, b を並べて a ≥ b と書く設計は採らない**——それは非対称性を消す。

## ★不定性を型の構造として主張しないこと(2026-08-14 の設計判断)

一度 `possibleThetaImages : Indeterminacy → Set Amb`(不定性で添字づけた族)という
設計を検討したが、**採らなかった**。原文は

> the holomorphic hull of the union of the possible images of a Θ-pilot object

と書くだけで、可能な像が**何かで径数づけられる**とは言っていない。族として書くと
原文にない構造を主張することになる。また (Ind1)(Ind2)(Ind3) は3つとも性質が違い
——(Ind1)(Ind2) は「〜が誘導する不定性」(作用)だが (Ind3) は
「upper semi-compatible」という**両立性の主張**であって作用ではない——
3つをどう合成するかを原文は書いていない。添字型を置けばその合成を我々が
決めることになる。`Set (Set Amb)` で受ければ**その問い自体が生じない**。

(Ind1)(Ind2)(Ind3) は「**なぜ像が複数ありうるか**」の理由として
`possibleThetaImages` の docstring に逐語で置く。原文が不定性を置いている位置
(Θ 側は subject、q 側は not subject)と一致する。

## 「the situation of Theorem 3.11」を列挙していないこと

原文 (IUTchIII p.173):
> Corollary 3.12. (Log-volume Estimates for Θ-Pilot Objects) Suppose
> that we are in the situation of Theorem 3.11.

仮説はこの一文だけで、**何が「the situation」に含まれるかは列挙されていない**
(Theorem 3.11 は物理 p.153–159 の7ページ)。ここで受けるのは
**Corollary 3.12 が名指しで読み出す分だけ**であり、`the situation` の全体ではない。
すなわち本 `structure` は「Theorem 3.11 の状況」の**下界**である。
-/

namespace ABC3.Interface.IUTchIII

open ABC3.Meta

/-- Corollary 3.12 が読み出すデータ。**各フィールドは原文の逐語に対応する**(G1)。

像を `Set Amb`(部分集合)としたのは、原文の証明(物理 p.174)が
`n,◦U_{j,vQ} ⊆ n,◦U^Q_{j,vQ}` と部分集合を扱い、
「the **union** of the possible images」が集合の和として書かれているため。 -/
structure PilotObjectData where
  /-- 多輻的表現の台。パイロット対象の像が住む場所。

  原文 (IUTchIII p.174):
  > in the multiradial representation of Theorem
  > 3.11, (i), which we do not regard as subject to the indeterminacies (Ind1), (Ind2),
  -/
  Amb : Type
  /-- `Amb` の位相。コンパクト性を述べるために要る。 -/
  topology : TopologicalSpace Amb
  /-- **対数体積の本体**。Proposition 3.9, (i) は対数体積を `𝔐(−) → ℝ` として定める
  ——`𝔐(−)` は**コンパクト**集合の集まりであり、値域は `ℝ` である。

  原文 (IUTchIII p.115):
  > — where we write “M[frak](−)” for the set of nonempty compact open subsets of

  (アルキメデス側は「compact closures of nonempty open subsets」。値域はどちらも `ℝ`。)

  ★2026-08-14: 以前は `logVol_ne_top_of_isCompact` という**posit**でこれを表していたが、
  それは仮定ではなく**定義域と値域**である。ここで型に写したので、
  「コンパクトなら `+∞` でない」は `Skeleton` 側の**定理**になった
  (`Skeleton.IUTchIII.logVol_ne_top_of_isCompact`)。 -/
  logVolCompact : ∀ U : Set Amb, @IsCompact Amb topology U → ℝ
  /-- **手続き正規化モノ解析的対数体積**。Θ 側・q 側の**両方**をこれで測る
  (同じ多輻的表現の中で比較するのだから、`logVol` は1つでなければならない)。

  値域が `WithTop ℝ`(= ℝ ⋃ {+∞})なのは原文どおり——Θ 側は `+∞` を排除していない:

  原文 (IUTchIII p.173):
  > − |log(Θ[ul2])| ∈ R[bb]
  -/
  logVol : Set Amb → WithTop ℝ
  /-- **正則包**。Θ 側にだけ現れる操作。

  原文 (IUTchIII p.173):
  > of the holomorphic hull [cf. Remark 3.9.5, (i)]
  > of the union of the possible images of a Θ-pilot object [cf. Definition 3.8, (i)],
  -/
  hull : Set Amb → Set Amb
  /-- **Θ-パイロット対象の「可能な像」の集まり**。

  原文 (IUTchIII p.173):
  > of the union of the possible images of a Θ-pilot object [cf. Definition 3.8, (i)],

  ★**なぜ「可能な像」が複数ありうるのか** — 不定性を被るから。それが原文の言い方:

  原文 (IUTchIII p.174):
  > multiradial representation of Theorem 3.11, (i), which we regard as subject to
  > the indeterminacies (Ind1), (Ind2), (Ind3) described in Theorem 3.11, (i), (ii).

  内訳(Theorem 3.11 での導入。ここでは展開せず、逐語だけ記録する):

  原文 (IUTchIII p.154):
  > (Ind1) the indeterminacies induced by the automorphisms of the procession
  > of D[scr]^-prime-strips Prc(n,◦D[frak]^_T);

  原文 (IUTchIII p.154):
  > (Ind2) for each v_Q ∈ V[bb]^non_Q (respectively, v_Q ∈ V[bb]^arc_Q), the indeterminacies induced
  > by the action of independent copies of Ism [cf. Proposition 1.2, (vi)]

  原文 (IUTchIII p.156):
  > (Ind3) as one varies m ∈ Z[bb], the isomorphisms of (a) are “upper semi-
  > compatible”, relative to the log-links of the n-th column of the LGP-
  -/
  possibleThetaImages : Set (Set Amb)
  /-- **q-パイロット対象の像**——単数。ここに不定性が現れないことが、
  原文の「which we do **not** regard as subject to」を型で表している。

  原文 (IUTchIII p.174):
  > for the procession-normalized mono-analytic log-volume of the image of a
  > q-pilot object [cf. Definition 3.8, (i)], relative to the relevant Kummer isomor-
  > phisms [cf. Theorem 3.11, (ii)], in the multiradial representation of Theorem
  -/
  qImage : Set Amb
  /-- `|log(q[ul2])|`。原文が実数として扱う量。 -/
  qAbs : ℝ
  /-- 原文 (IUTchIII p.174):
  > In particular, |log(q[ul2])| > 0 is easily computed
  > in terms of the various q-parameters of the elliptic curve E_F
  -/
  qAbs_pos : 0 < qAbs
  /-- q 側の対数体積は `−|log(q[ul2])|` であり、しかも**実数**(`+∞` ではない)。
  Θ 側の値域(`ℝ ⋃ {+∞}`)との違いは原文由来である。

  原文 (IUTchIII p.174):
  > Write
  > − |log(q[ul2])| ∈ R[bb]
  > for the procession-normalized mono-analytic log-volume of the image of a
  > q-pilot object [cf. Definition 3.8, (i)],
  -/
  qLogVol_eq : logVol qImage = ((-qAbs : ℝ) : WithTop ℝ)
  /-- **可能な出力対数体積の集まり**(原文の `ℝ_{≤−|log(Θ̲̲)|}`)。

  ★2026-08-14 追加。反証(`cor_3_12_refutable_under_current_interface`)を殺すために
  原文へ当たり直して見つけたもの。**Corollary 3.12 の statement にはなく、
  その証明の Step (xi-e) にある**(物理 p.184)。

  原文 (IUTchIII p.184):
  > The multiradial construction algorithm of Theorem 3.11, followed by for-
  > mation of the holomorphic hull and application of the log-volume, yields a
  > collection of possible log-volumes of pilot-object output data
  -/
  outputLogVolumes : Set ℝ
  /-- ★2026-08-15: ここにあった **原典の主張** は `Skeleton/IUTchIII/Cor312Claims.lean` へ出した。

  出した理由: 主張を `Interface` のフィールドにすると**仮説**になるので、
  依存グラフの辺の先が**展開不能**になる(掘るべき対象が goal ではなく hypothesis になる)。
  出したのは次の5つ:
  `PossibleImagesContained`(Ob1, p.131) / `LogShellPacketCompact`(p.153, p.31, p.146) /
  `HullCompactOfRelCompact`(Remark 3.9.5, (i), p.127) /
  `OutputLogVolumesEq`(Step (xi-e), p.184) / `QLogVolMem`(Step (xi-f), p.184)。

  ★出した先は `theorem ... := sorry` **ではなく** 名前付きの `Prop` である。
  `PilotObjectData` はデータの袋であって原典の対象を同定する条件を持たないので、
  これらは任意の `D` については**偽**であり、`sorry` で置けば嘘になる。
  反証は `Check/IUTchIII/Cor312Degenerate.lean` で実際に構成した。

  以下このフィールドは、その事実の記録としてのみ置く(データではない)。 -/
  claimsMovedToSkeleton : True
  /-- **log-shell(モノ解析的整構造)** `I(…)`。Theorem 3.11, (i), (a) が
  多輻的表現のデータとして挙げるもの。

  原文 (IUTchIII p.153):
  > analytic integral structures

  ## ★このフィールド群の由来 — 1段深い転写(2026-08-14)

  一度ここには **`thetaUnion_isCompact`**(可能な像の和集合はコンパクト)という
  posit が1本あった。出所は物理 p.175:

  原文 (IUTchIII p.175):
  > pactness of the 1,◦U_j,v_Q [where j ∈ |F_l|, v_Q ∈ V[bb]_Q], together with the definition of the
  > log-volume, that the quantity − |log(Θ[ul2])| is finite, hence negative

  ★★**バーの有無に注意**(400 dpi 目視): 原文がコンパクトだと言っているのは
  `¹,°𝒰_{j,v_ℚ}`(**オーバーバー無し**)であり、p.174 の定義により
  これは「the various **unions** … of the **possible images**」——すなわち
  **可能な像の和集合**である。正則包の方は `¹,°𝒰̄`(バー有り)で別物。
  `pdftotext` はオーバーバーを落とすので `.txt` では区別できない。

  ★原文は「[easily verified]」と書くだけで**検証を書いていない**(実測:
  論文全体で `compact` は 33 件、`¹,°𝒰` のコンパクト性を確立する箇所は 0 件)。

  そこで1段深く分解した——`thetaUnion_isCompact` を posit するのをやめ、
  下の4つ(`logShellPacket` / `logShellPacket_isCompact` /
  `possibleThetaImages_subset_logShellPacket` / `hull_isCompact_of_subset_isCompact`)と
  `logVol_ne_top_of_isCompact` に置き換え、有限性は `Skeleton` で**導出**する。
  測定結果は `Check/IUTchIII/Cor312Degenerate.lean`。 -/
  logShellPacket : Set Amb
  /-- `logVol` はコンパクト領域上で `logVolCompact` と一致する
  ——すなわち `logVol` は Proposition 3.9, (i) の対数体積の**拡張**である。

  Corollary 3.12 が `−|log(Θ̲̲)| ∈ ℝ ⋃ {+∞}` と書く以上、拡張は要る。
  拡張が非コンパクト側でどう振る舞うかは、`Interface` では**指定しない**
  (Remark 3.9.5, (i) は「相対コンパクトでなければ正則包は `I^ℚ` 全体」と述べ、
  Remark 3.9.7, (ii) は相対コンパクト側の拡張が `−∞` を許すと述べる。
  我々はどちらも写していない)。 -/
  logVol_eq_of_isCompact : ∀ (U : Set Amb) (h : @IsCompact Amb topology U),
    logVol U = ((logVolCompact U h : ℝ) : WithTop ℝ)

/-- Track B は何を作らねばならないか。

★**退化 witness なら今すぐ作れる**(`Check/IUTchIII/Cor312Degenerate.lean` で実際に作った)。
にもかかわらず `nonvacuous` ではなく `waiting` を置くのは、退化 witness で G2 を
満たしても作業キューから消えるだけで何も進まないから——`tools/check.mjs` 冒頭 B5 の穴。 -/
def PilotObjectData.waiting : WaitingFor :=
  { what := "Θ-パイロット対象 / q-パイロット対象(Definition 3.8, (i))・正則包(Remark 3.9.5, (i))・手続き正規化モノ解析的対数体積(Proposition 3.9)・不定性 (Ind1)(Ind2)(Ind3)(Theorem 3.11)"
    trackB := "Found/IUTchIII — その手前に mono-theta 環境・log-shell・Frobenioid・tempered 基本群が要る(mathlib に 0 件、PLAN 事実3)" }

end ABC3.Interface.IUTchIII
