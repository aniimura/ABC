import ABC3.Skeleton.PGC.Section2
import ABC3.Skeleton.PGC.Section3Defs

/-!
# [pGC] §3 — 命題

設定・記号は `ABC3/Skeleton/PGC/Setup.lean`。§2 の `RamificationFiltration`・
`RecoverableAsAddModule` を継続して使う。

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic
Local Fields/section-3.html`(PDF 目視確認 2026-09-03、物理 p.6)。

## ★2026-09-08: 定義を `Section3Defs.lean` へ割った

`filteredGroupOf` と `IsUniformizing`(どちらも**定義**)は
`ABC3/Skeleton/PGC/Section3Defs.lean` へ移した。本ファイルには**定理だけ**が残る
(`Section1Defs.lean`・`Section2Defs.lean` と同じ作法)。
★名前空間は同じ `ABC3.Skeleton.PGC` なので完全修飾名は 1 文字も変わらない。
理由(`Found/PGC/` の 200 本が本ファイルを import できなかったこと)は
`Section3Defs.lean` の docstring に測定つきで書いた。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Corollary 3.1 -/

/-- **[pGC] Corollary 3.1**

原文 (pGC p.6):
> Given a continuous Q_p[Γ_K]-vector space of finite Q_p-dimension, the issue of whether or
> not V is Hodge-Tate (as well as the invariants d_V(i)) can be determined entirely
> group-theoretically from the filtered group Γ_K.

## 形式化上の簡略化(逸脱として記録)

「Hodge-Tate である」の定義(`d_V ≤ dim_{Q_p}(V)`、`d_V(i)` は `V(−i)⊗K̄^` の Γ_K-不変部分の
次元)は Proposition 2.2 の `CompKbar`(K̄^ の抽象化)とテンソル積・Tate twist を要し、
本セッションでは未構築。ここでは「Hodge-Tate である」を**抽象的な述語パラメータ**
`isHodgeTate : ∀ K, Prop` として受け取り、その値が α で保たれることだけを主張する
——実際の Hodge-Tate 述語の構成は独立した課題として残す。

## 逸脱の訂正(2026-09-04)

★原文は明示的に「the **filtered** group Γ_K」(フィルトレーション込みの群)から
回復できると述べている。以前の形式化は `_α` を**裸の**
`ContinuousMulEquiv K.absGal K'.absGal`(フィルトレーションと無関係な同型)に
取っていた——`_RF` はパラメータとして受け取りながら `_α` の型を一切制約していない、
という不整合(先頭の `_` は「未使用」の意味そのままだった)。

数論的に見ても、これは看過できない差である: `Γ_K` を**裸の**副有限群として見た
同型類は(奇素数 p では)`p` と `[K:Q_p]` だけで決まり、`K` 自身には依らない
(Iwasawa の型の古典的事実)——したがって裸の `α` が存在するだけでは
`K ↦ K` に依存する述語(`isHodgeTate` のような自由なパラメータ)の不変性を
導く根拠に**なりえない**。原文がわざわざ「filtered group」と言っているのは
まさにこの理由による。

**訂正**: `_α` の型を `FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K')`
に直した(`Theorem 4.2` が `FilteredGroup.OuterIso` を使うのと同じ発想——
高次分岐群のフィルトレーションを保つ同型のみを許す)。`RF` は今回から実際に
型の中で使われる(`_RF` → `RF`)。

## 依拠する境界外の結果

- [1] Serre, *Abelian ℓ-adic Representations*, Chapter III §1.2(d_V ≤ dim_{Q_p}(V)、
  Hodge-Tate の定義そのもの)。mathlib 不在。
- Proposition 2.2(K̄^ の回復)への直接依存。

## ★★★2026-09-05: `V` が**使われていなかった**(訂正)

旧形は `V : PAdicLocalField p → Type*` を仮説に取りながら、結論
`isHodgeTate K ↔ isHodgeTate K'` に **`V` が一度も現れない**——
すなわち「`K` 添字の**任意の** Prop 族は filtered iso で不変」と
言っていた。これは `Check/PGC/Theorem42Degenerate.lean`(自由な `Φ`)・
`Check/PGC/Cor33Degenerate.lean`(無関係な `ρ`・`ρ'`)と同じ型の欠陥である
——ただし反例には `K ≠ K'` の witness が要るので、今のところ**反証は
できていない**(`Check/PGC/RefutationAttempts.lean` の壁)。

**修理**: 原文は「**与えられた** `V` が Hodge-Tate かどうか」を言うのだから、
`isHodgeTate` を **`V` にも依存する述語**にし、`α` で対応する2つの作用を
持つ**同一の** `V` について比べる形にした(`Cor 3.3` の修理と同じ発想:
移送で結ばれた対を比べる)。`K' = K`・`α = 恒等` では `_hcompat` が
`sK' = sK` を強制するので、`Iff.rfl` に落ちる。

## ★★検討中の懸念(2026-09-04、訂正込み——上の修理で `V` は使われるようになった)

`isHodgeTate : ∀ K, Prop` は**完全に自由なパラメータ**——`V` の構造から
導かれる制約が一切無い。一見 `Cor 1.3` の `I_K`・`Prop 1.2` の `RD` の
「自由なデータによる退化」と同じ型の罠に見えるが、精査すると**事情が違う**:
`K′=K` の場合(`_α` がどんな `FilteredGroup.Iso Γ_K≅Γ_K` であっても)は
`isHodgeTate K ↔ isHodgeTate K` が `Iff.rfl` で**自明に成り立つ**——
`isHodgeTate` が Galois 群の元(共役の選び方)に一切依存しないパラメータ
だから、`Prop 1.1` の共役不変性(`cyclotomicCharacter_conj_invariant`)の
ような別途の証明さえ要らない。ゆえに本当の障害は `K≠K′` かつ
`FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K')` が
存在するという**witness の構成そのもの**に尽きる——これは
`Check/PGC/RefutationAttempts.lean` が pGC §1 の3定理について
既に詳細に調べ尽くした障害(「異なる2つの p進局所体の間の連続同型」が
現状の道具では構成できない、0/3 反証できなかった)と**同一の根本原因**
であり、`isHodgeTate` 固有の新しい罠ではない。したがって本項目は
`Prop 1.1`・`Prop 1.2`・`Cor 1.3` と同じ意味で「反証もできないし証明も
できない」——未解決のまま、`isHodgeTate` を制約する訂正は不要と見る。 -/
theorem cor_3_1 (RF : RamificationFiltration p)
    (isHodgeTate : ∀ (K : PAdicLocalField p) (V : Type)
      [AddCommGroup V] [Module ℚ_[p] V] [SMul K.absGal V], Prop) :
    ∀ {K K' : PAdicLocalField p}
      (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
      (V : Type) [AddCommGroup V] [Module ℚ_[p] V] [FiniteDimensional ℚ_[p] V]
      (sK : SMul K.absGal V) (sK' : SMul K'.absGal V)
      (_hcompat : ∀ (g : K.absGal) (x : V), sK'.smul (α.equiv g) x = sK.smul g x),
      @isHodgeTate K V _ _ sK ↔ @isHodgeTate K' V _ _ sK' := sorry

def cor_3_1.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def cor_3_1.needs : List ProofObligation :=
  [ .citation "[1] Serre, Abelian ℓ-adic Representations, Chapter III §1.2"
      "Hodge-Tate の定義(d_V ≤ dim_{Q_p}(V))"
      (.absent "mathlib v4.31.0-rc2 実測: Hodge-Tate 表現論に相当する宣言はゼロ件。★2026-09-06 に再測: re:`HodgeTate|hodgeTate|Hodge-Tate|isHodgeTate|hodgeTateWeight`→0") 6,
    .otherPaper "pGC" "Proposition 2.2" 5 ]

/-! ## Corollary 3.3 -/

/-- **[pGC] Corollary 3.3**

原文 (pGC p.6):
> Given a continuous E[Γ_K]-module V of E-dimension 1, the issue of whether or not V is
> uniformizing can be determined entirely group-theoretically from the filtered group Γ_K.

## 逸脱の訂正(2026-09-04)

Corollary 3.1 と同じ理由(原文「the filtered group Γ_K」、裸の同型では
`K` に依存する述語の不変性を導く根拠にならない——`Γ_K` の裸の同型類は
`p` と `[K:Q_p]` だけで決まるという古典的事実)で、`_α` の型を
`FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K')` に訂正した
(`RF` も `_RF` → `RF` に)。

## 依拠する境界外の結果

- [1] Serre, Chapter III, Appendix §5(d_V(i) による uniformizing の判定条件式)。
- Corollary 3.1 への直接依存(d_V(i) が回復できることを使う)。

## ★★★2026-09-05: `ρ` と `ρ'` が無関係なのは**やはり罠だった**(訂正)

下の 2026-09-04 の見立て(「反証はやはり効かない」)は**誤りだった**。
`K' = K`・`α = 恒等` を取ると、旧形は「どんな2つの表現も uniformizing 性が
一致する」と言っていることになる。そして一致しない2つは実際に作れる:

* **真になる側**: `E := K`・`ι := id`・`I := U_K` とし、
  `Γ_K ↠ 𝒪_K^×`(`Found/PGC/AbsGalUnitsSurjective.lean`、Lubin-Tate 相互律を
  無条件化したもの)から `ρ := Γ_K ↠ 𝒪_K^× ↪ K^×`、`toGal := その切断`。
* **偽になる側**: 自明な表現 `ρ' = 1`(`Check/PGC/Def32Degenerate.lean::
  not_isUniformizing_one`)。

証明は `Check/PGC/Cor33Degenerate.lean::cor_3_3_statement_false`(`sorry` 無し)。
2026-09-04 の見立てが外れたのは、当時「真になる `ρ` を作るには相互律相当が
要る」ところで止まっていたからで、その相互律の半分が構築された今は作れる。

**修理**: 原文は**ひとつの** `V` を `α` で移して見比べるのだから、
`ρ'` を `α` で `ρ` と結ぶ条件 `∀ g, ρ' (α.equiv g) = ρ g` を課した。
`K' = K`・`α = 恒等` では `ρ' = ρ` が強制されるので反例は塞がる。

## ★★検討中の懸念(2026-09-04、訂正込み——上のとおり**この見立ては外れた**)

`ρ : K.absGal →* Eˣ`・`ρ' : K'.absGal →* Eˣ` が `_α` と無関係に自由に
選べる点は一見 `Cor 3.1` と同じ罠に見えるが、精査すると**反証はやはり
効かない**——`toGal` が固定されていても、`IsUniformizing K E toGal ρ`
自体が `ρ`(と `toGal`)の**具体的な組**に強く依存する非自明な存在命題
(開集合 `I` 上で `ρ∘toGal` がある体準同型 `ι` と一致する)であり、
`ι` は単射(体準同型)ゆえ `I` 上で単射的に変化しなければならない。
`ρ∘toGal` がその条件を満たす——すなわち `IsUniformizing` が**真になる**
——ような `(toGal,ρ)` の組を構成するには、それ自体が局所類体論の
相互律に相当する非自明な数学的内容を要る(`toGal` を「悪く」選べば
両辺とも恒等的に偽になり ↔ は自明に真、「良い」`toGal` を選ぶには
相互律相当の構成が要る、という板挟み)。ゆえに `Cor 3.1` と同じ根本
原因(`Check/PGC/RefutationAttempts.lean` が突き止めた、局所類体論を
経由しない反証・証明のどちらも現状の道具では届かないという壁)に帰着する
——`ρ,ρ'` の非拘束それ自体は独立した罠ではない。 -/
theorem cor_3_3 (RF : RamificationFiltration p)
    (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : ∀ K : PAdicLocalField p, {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) :
    ∀ {K K' : PAdicLocalField p}
      (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
      (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ)
      (_hρ : ∀ g : K.absGal, ρ' (α.equiv g) = ρ g),
      IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ' := sorry

def cor_3_3.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.3", sectionId := "cor-3-3" }

def cor_3_3.needs : List ProofObligation :=
  [ .citation "[1] Serre, Abelian ℓ-adic Representations, Chapter III, Appendix §5"
      "d_V(1) = [E:K]; d_V(0) = [E:K]·([K:Q_p]−1) による uniformizing の判定"
      (.absent "mathlib v4.31.0-rc2 実測: 該当なし(Hodge-Tate 重み d_V(0)・d_V(1) による uniformizing の判定)。★2026-09-06 に再測: re:`HodgeTate|hodgeTate|Hodge-Tate|isHodgeTate|hodgeTateWeight`→0") 6,
    .otherPaper "pGC" "Corollary 3.1" 6 ]

end ABC3.Skeleton.PGC
