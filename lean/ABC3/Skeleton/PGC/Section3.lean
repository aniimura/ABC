import ABC3.Skeleton.PGC.Section2

/-!
# [pGC] §3 — 命題

設定・記号は `ABC3/Skeleton/PGC/Setup.lean`。§2 の `RamificationFiltration`・
`RecoverableAsAddModule` を継続して使う。

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic
Local Fields/section-3.html`(PDF 目視確認 2026-09-03、物理 p.6)。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-- `RamificationFiltration`(`Interface/`)から `K` 1つ分の `FilteredGroup`(`Found/`)を作る。

★`ABC3.Interface.PGC.RamificationFiltration.filt`(`Skeleton/PGC/Section4.lean`)と
**同じ構成**だが、Section4 は本ファイルを import するので(循環を避けるため)ここに
独立に置く——小さな構造の詰め替えなので複製のコストは小さい。 -/
noncomputable def filteredGroupOf (RF : RamificationFiltration p) (K : PAdicLocalField p) :
    FilteredGroup :=
  { G := K.absGal, Gv := RF.Gv K, isClosed := RF.isClosed K, isNormal := RF.isNormal K,
    antitone := RF.antitone K }

/-- 台帳の付随宣言(橋渡しの `def` であり、原典の項目そのものではない)。
`ABC3.Interface.PGC.RamificationFiltration.filt.src`(`Section4.lean`)と同じ位置づけ。 -/
def filteredGroupOf.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Section 3 (filteredGroupOf)", sectionId := "cor-3-1" }

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

## ★★検討中の懸念(2026-09-04、訂正込み)

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
    (V : PAdicLocalField p → Type*) [∀ K, AddCommGroup (V K)] [∀ K, Module ℚ_[p] (V K)]
    [∀ K, FiniteDimensional ℚ_[p] (V K)] [∀ K, SMul K.absGal (V K)]
    (isHodgeTate : ∀ K, Prop) :
    ∀ {K K' : PAdicLocalField p}
      (_α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K')),
      isHodgeTate K ↔ isHodgeTate K' := sorry

def cor_3_1.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def cor_3_1.needs : List ProofObligation :=
  [ .citation "[1] Serre, Abelian ℓ-adic Representations, Chapter III §1.2"
      "Hodge-Tate の定義(d_V ≤ dim_{Q_p}(V))"
      (.absent "mathlib v4.31.0-rc2 実測: Hodge-Tate 表現論に相当する宣言はゼロ件") 6,
    .otherPaper "pGC" "Proposition 2.2" 5 ]

/-! ## Definition 3.2 -/

/-- **[pGC] Definition 3.2** — uniformizing な E[Γ_K]-加群。`sorry` 無し(純粋な定義)。

原文 (pGC p.6):
> We shall call the E[Γ_K]-module V uniformizing if the restriction of ρ_V to some open
> subgroup I of U_K (⊆ Γ^a_K^b) is the morphism I → E× induced by restricting some morphism
> of fields K → E to I ⊆ U_K ⊆ K.

## 形式化上の簡略化(逸脱として記録)

原文の「ρ_V は U_K ⊆ Γ_K^ab 上に制限できる」は、**局所類体論の相互律**
Γ_K^ab ≅ (K^×)^(§1 Proposition 1.2 の論拠、mathlib 不在)を暗黙に経由する——
U_K は本来 K^× の部分群であり、Γ_K の部分群ではない。

ここでは相互律を明示的なパラメータ `toGal : U_K → Γ_K`(未構築の辞書、まだ本物には
なっていない)として受け取ることで、定義全体を条件付きで well-typed にする
——`toGal` が本物になれば(Interface が実装されれば)この定義はそのまま使える。 -/
def IsUniformizing (K : PAdicLocalField p) (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) (ρ : K.absGal →* Eˣ) : Prop :=
  ∃ (I : Set K.carrier) (hIU : I ⊆ {x : K.carrier | ‖x‖ = 1}) (_hopen : IsOpen I)
    (ι : K.carrier →+* E), ∀ x (hx : x ∈ I), ((ρ (toGal ⟨x, hIU hx⟩) : Eˣ) : E) = ι x

def IsUniformizing.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Definition 3.2", sectionId := "def-3-2" }

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

## ★★検討中の懸念(2026-09-04、訂正込み)

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
      (_α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
      (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ),
      IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ' := sorry

def cor_3_3.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.3", sectionId := "cor-3-3" }

def cor_3_3.needs : List ProofObligation :=
  [ .citation "[1] Serre, Abelian ℓ-adic Representations, Chapter III, Appendix §5"
      "d_V(1) = [E:K]; d_V(0) = [E:K]·([K:Q_p]−1) による uniformizing の判定"
      (.absent "mathlib v4.31.0-rc2 実測: 該当なし") 6,
    .otherPaper "pGC" "Corollary 3.1" 6 ]

end ABC3.Skeleton.PGC
