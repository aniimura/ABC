import ABC3.Interface.IUTchIII.PilotObjects
import ABC3.Meta.Claim

/-!
# Corollary 3.12 が依拠する **原典の主張**(2026-08-15 に `Interface` から出した)

## なぜ出したか

これらは 2026-08-14 まで `Interface` の `PilotObjectData` の**フィールド**だった。
フィールドは仮説なので、依存グラフの**辺の先が展開不能**になる——
掘るべき対象が goal ではなく hypothesis になり、足場が消える。
実際、`tools/check.mjs` の依存グラフ印字で
`[IUTchIII] Theorem 3.11` / `Remark 3.9.5` が **未展開の葉**のまま動かなかった。

## ★なぜ `theorem ... := sorry` にしなかったか(設計判断)

「主張は `Skeleton` に `sorry` で置く」という素朴な形は**採れない**。
`PilotObjectData` は**データの袋**であって、それが原典の対象であることを
同定する条件を持っていない。したがって

```
theorem ob1 (D : PilotObjectData) : ∀ U ∈ D.possibleThetaImages, U ⊆ D.logShellPacket
```

は**任意の `D` については偽**である。`sorry` で置けば、それは
「いつか埋まる借金」ではなく**永久に埋まらない偽の言明**になる。
posit(仮説)は「仮定する」と正直に言っているが、偽の `sorry` は嘘である。

★これは意見ではない。`Check/IUTchIII/Cor312Degenerate.lean` の
`ob1_is_refutable` / `packetCompact_is_refutable` で**反証を構成して確かめた**。

そこで採った形は: **主張に名前を付けた `Prop` を `Skeleton` に置き、
`cor_3_12` は明示的な仮説として受け取る**。これで

- 辺の先が展開可能になる(`.src` と `.needs` を持つ**節点**になる)
- 依存が `.needs` として見える(中層はここに現れる)
- **偽の言明を主張しない**
- 借金が `cor_3_12` の**statement に見える**(フィールドは隠れていた)

★★**`sorry` は増えていない。** それは減らそうとしたからではなく、
偽を `sorry` で置くことを拒んだからである。「`sorry` が増えるのが正しい」という
見立ては、**同定条件を書ける場合に限って**正しい。
-/

namespace ABC3.Skeleton.IUTchIII

open ABC3.Interface.IUTchIII ABC3.Meta

/-! ### 有限性を支える3つの主張 -/

/-- **(Ob1)** 可能な像は log-shell のテンソルパケットに含まれる。

原文 (IUTchIII p.131):
> various “possible images” that occur as the output of the multiradial al-
> gorithms under consideration are regions — i.e., in essence, elements ∈P
> — contained in tensor packets of log-shells I[scr]_k

★2026-08-14 訂正の記録: 一度これを「原文に無い」と報告したが**誤り**だった。
Lean 側で付けた名前の語(`subset` / `contains`)で grep していたためで、
原文は原文側の語(`possible images`)を主語にしていた。 -/
def PossibleImagesContained (D : PilotObjectData) : Prop :=
  ∀ U ∈ D.possibleThetaImages, U ⊆ D.logShellPacket

def PossibleImagesContained.src : Source :=
  { paper := "IUTchIII", pdfPage := 131, item := "Remark 3.9.5, (vii), (Ob1)",
    sectionId := "remark-3-9-5-vii-ob1" }

def PossibleImagesContained.needs : List ProofObligation :=
  [ .implicitStep
      "★2026-08-15 に辺を1本**落とした**: 一度 `[IUTchIII] Theorem 3.11`(物理 p.153)への辺を書いていたが、原文 p.131 は「arises from applying the multiradial algorithms of Theorem 3.11 **below**」と書いており、**前方参照**である。Remark 3.9.5(物理 p.126-131)は Theorem 3.11(物理 p.153)より**前**にあるので、その主張が Theorem 3.11 に依拠することはできない。★迷いを記録する: Ob1 の内容は多輻的アルゴリズムの出力について述べているので、内容としては Theorem 3.11 を前提しているとも読める。それでも辺にしないのは、(a) 原文が `below` と明示していること、(b) `cor_3_12` から Theorem 3.11 への辺が既にあり依存は表現されていること、の2点による" 131,
    .implicitStep
      "log-shell の**テンソルパケット** `ℐ_k` の構成そのもの(j 上のテンソル積)を写していない。我々は `logShellPacket : Set Amb` という単一の部分集合で受けている" 131,
    .implicitStep
      "原文の `elements ∈ P`(領域のなす集合 `P`)を写していない。我々は `Set (Set Amb)` で受けており、領域の圏構造は落としている" 131 ]

/-- **log-shell のテンソルパケットはコンパクト**。

原文 (IUTchIII p.153):
> analytic integral structures

パケットが Theorem 3.11, (i), (a) の多輻的表現のデータとして入る。
コンパクト性そのものは別の箇所で述べられる:

原文 (IUTchIII p.31):
> satisfies the following properties: (anon) I†Fv is compact, hence of finite log-

原文 (IUTchIII p.146):
> O(−) = I((−)) ⊆IQ((−))
-/
def LogShellPacketCompact (D : PilotObjectData) : Prop :=
  @IsCompact D.Amb D.topology D.logShellPacket

def LogShellPacketCompact.src : Source :=
  { paper := "IUTchIII", pdfPage := 153, item := "Theorem 3.11, (i), (a)",
    sectionId := "thm-3-11-i-a" }

def LogShellPacketCompact.needs : List ProofObligation :=
  [ .otherPaper "[AbsTopIII]" "Corollary 5.10, (i)(log-shell の基本性質)" 145,
    .implicitStep
      "★2026-08-15 に辺を**差し替えた**: 一度 `[AbsTopIII] (L1)`(物理 p.5)へ向けていたが、原文 [IUTchIII] 物理 p.31 が実際に引用しているのは「is compact, hence of finite log-volume [cf. [AbsTopIII], Corollary 5.10, (i)]」である。(L1) は導入部の要約であって、それ自身 [cf. Corollary 5.10, (i)] と本体を指している。導入部から本体への指しは辺にしない(導入部は証明を持たない)" 31,
    .implicitStep
      "★単一の log-shell のコンパクト性(p.31)から**テンソルパケット**のコンパクト性への段を写していない。有限テンソル積であることが要るが、原文はここを明示していない" 146 ]

/-- **相対コンパクトな領域の正則包はコンパクト**。

原文 (IUTchIII p.127):
> If αU (respectively, AU; A,αU) is relatively compact, then we define
> the holomorphic hull of αU (respectively, AU; A,αU) to be the smallest subset of
-/
def HullCompactOfRelCompact (D : PilotObjectData) : Prop :=
  ∀ U K : Set D.Amb, @IsCompact D.Amb D.topology K → U ⊆ K →
    @IsCompact D.Amb D.topology (D.hull U)

def HullCompactOfRelCompact.src : Source :=
  { paper := "IUTchIII", pdfPage := 127, item := "Remark 3.9.5, (i)",
    sectionId := "remark-3-9-5-i-hull" }

def HullCompactOfRelCompact.needs : List ProofObligation :=
  [ .implicitStep
      "原文は正則包を「the smallest subset of …」と**最小性**で定義する。その最小集合の**存在**を我々は写していない(`hull : Set Amb → Set Amb` という演算として受けている)" 127,
    .implicitStep
      "相対コンパクトでない枝(正則包が `I^ℚ` 全体になり対数体積が `+∞` になる場合)を写していない。Corollary 3.12 の有限性はコンパクトな枝の側だけに乗っている" 127 ]

/-! ### 不等式を支える2つの主張(証明 Step (xi-e)(xi-f)) -/

/-- **(xi-e)** 可能な出力対数体積の集まりは「Θ 側の対数体積以下の実数」全体である。

原文 (IUTchIII p.184):
> The multiradial construction algorithm of Theorem 3.11, followed by for-
> mation of the holomorphic hull and application of the log-volume, yields a
> collection of possible log-volumes of pilot-object output data

原文 (IUTchIII p.184):
> R[bb]_≤−|log(Θ[ul2])| ⊆ R[bb]
-/
def OutputLogVolumesEq (D : PilotObjectData) : Prop :=
  D.outputLogVolumes
    = {x : ℝ | (x : WithTop ℝ) ≤ D.logVol (D.hull (⋃₀ D.possibleThetaImages))}

def OutputLogVolumesEq.src : Source :=
  { paper := "IUTchIII", pdfPage := 184,
    item := "Corollary 3.12 (proof, Step (xi-e), (xi-f))",
    sectionId := "cor-3-12-step-xi" }

def OutputLogVolumesEq.needs : List ProofObligation :=
  [ .otherPaper "[IUTchIII]" "Theorem 3.11(多輻的構成アルゴリズム)" 153,
    .implicitStep
      "原文が `⊆ ℝ` と書いていることは、`−|log(Θ)|` 自身が実数であることを含意しない(`−|log(Θ)| = +∞` でも `ℝ_{≤+∞} = ℝ ⊆ ℝ`)。有限性は別に要る" 184 ]

/-- **(xi-f)** q 側の対数体積はその集まりに属する。原文はここから不等式が
「then follows formally」だと述べる。

原文 (IUTchIII p.184):
> The inclusion − |log(q[ul2])| ∈ R[bb]_≤−|log(Θ[ul2])|, hence also the inequality

★★**落とした限定(意図的。記録)**: 直前の文はこの帰属を条件つきで、しかも濁して述べる。

原文 (IUTchIII p.184):
> log-volumes of output data is subject to the condition that this con-
> struction of output data possibilities constitutes, in particular, a construc-
> tion [perhaps only up to some sort of “approximation”, as a result of vari-
> ous indeterminacies] of the pilot-object log-volume of the input data

我々が写したのは**無条件の**帰属である。`approximation` は原文のどこにも
量化されていないため、写せば濁しの強さを我々が決めることになる。 -/
def QLogVolMem (D : PilotObjectData) : Prop := (-D.qAbs) ∈ D.outputLogVolumes

def QLogVolMem.src : Source :=
  { paper := "IUTchIII", pdfPage := 184,
    item := "Corollary 3.12 (proof, Step (xi-f))",
    sectionId := "cor-3-12-step-xi" }

def QLogVolMem.needs : List ProofObligation :=
  [ .implicitStep
      "原文 p.184 の限定『subject to the condition』と『perhaps only up to some sort of approximation, as a result of various indeterminacies』を写していない。仮説の強化にあたりうる意図的な単純化(check.mjs 冒頭 A6)" 184,
    .implicitStep
      "★2026-08-15 に辺を1本**落とした**: 一度 `[IUTchII] Corollary 4.10, (i)`(物理 p.158)への辺を書いていたが、**原文 p.184 に `[IUTchII]` の引用は 0 件**である(実測)。あれは私が『出力データの構成が入力データの pilot-object 対数体積を与える段が依拠するはずだ』と**推測して書いたもの**であり、原文が言っていない依存だった。なお `cor_3_12` から `[IUTchII] Corollary 4.10, (i)` への辺は原文 p.174(記号 △ の定義)にあるので、そちらは残る" 184 ]

end ABC3.Skeleton.IUTchIII
