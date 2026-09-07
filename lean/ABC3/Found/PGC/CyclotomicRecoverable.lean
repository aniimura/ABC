import ABC3.Found.PGC.CyclotomeTransport
import ABC3.Found.PGC.CyclotomicFromAbelianization
import ABC3.Found.PGC.AdjoinFieldClosure
import Mathlib.FieldTheory.Galois.Infinite
import Mathlib.Topology.Algebra.Group.TopologicalAbelianization

/-!
# 経路 Λ11 —— 群論的円分子 `Λ_n(H)` は `μ_{p^n}` である(同型として)

[pGC] Proposition 1.1 への**経路 Λ**の最終節点。Λ1
(`Found/PGC/TopAbelianization.lean`)・Λ2(`Found/PGC/TorsionCyclotome.lean`)・
Λ10(`Found/PGC/CyclotomeTransport.lean`)が用意した**群論側**の円分子と、
Λ8(`Found/PGC/ArtinMap.lean`)・Λ9(`Found/PGC/CyclotomicFromAbelianization.lean`)が
用意した**局所類体論側**の `tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)` を、
ここで初めて突き合わせる。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

## 何が埋まって、何が残ったか(★これが本ファイルの主眼)

Λ10 は Proposition 1.1 を `TorsionCyclotomeIsCyclotomic p`(1 本の補題)に還元していた。
その 1 本は 2 つの内容の合成である:

| 段 | 内容 | 本ファイルの結果 |
|---|---|---|
| (i) | `Λ_n(S) ≅ ℤ/p^n`(**同型があること**) | ★**証明した**(`nonempty_cyclotomeEquivZMod`) |
| (ii) | 共役作用が `χ_{F,n}` 倍(**作用の仕方**) | ★**残った**(`CyclotomeConjIsCyclotomic`) |

すなわち Proposition 1.1 は、いまや

  `CyclotomeConjIsCyclotomic p → (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal`
  (`cyclotomicCharacter_recoverable_of_cyclotomeConj`)

という形で、**存在量化子を 1 つも含まない 1 文**に還元されている:

  「`μ_{p^n} ⊆ L_S` なる開正規 `S ⊴ Γ_F` の円分子 `Λ_n(S)` の上で、
   `g ∈ Γ_F` の共役作用は `x ↦ x^{χ_{F,n}(g)}` である。」

★Λ10 の `TorsionCyclotomeIsCyclotomic` は「そのような同型 `e` が存在する」という
存在命題だったので、(i) と (ii) が混ざっていた。本ファイルは (i) を実際に構成して
取り除き、残りを (ii) だけにした。

## (i) の道筋(4 段、いずれも在庫の合成)

1. **抽象核**(`topAbelianizationEquivGalFixedField`、分岐・付値・局所体の語彙ゼロ):
   `E/k` が Galois なら、位相的アーベル化 `Gal(E/k) ⧸ [Gal,Gal]‾` は
   `Gal(E^{[G,G]‾} / k)` と同型。mathlib の
   `InfiniteGalois.normalAutEquivQuotient`(閉正規部分群についての無限次 Galois 対応)に
   `N := (commutator _).topologicalClosure` を代入するだけ。
2. **具体層**(`topAbelianizationEquivAbelianGal`): `k := K`、`E := K̄` を代入すると
   `abelianClosure K = fixedField ((commutator Γ_K).topologicalClosure)`(定義)なので
   `Γ_K^{ab}(位相的) ≅ Gal(K^ab/K)` を得る。★**これが Λ1/Λ2 と Λ8/Λ9 を繋ぐ 1 本**である。
3. **Λ9 を当てる**: `tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)`
   (`exists_abelianGalTorsion_equiv_rootsOfUnity`、仮説ゼロ)。
   `cyclotome A m` と `powTorsion (A^{ab}) m` は**定義が同じ**(どちらも
   `(powMonoidHom m).ker`)なので橋は要らない。
4. **開部分群を固定体に移す**: `absGalFixedFieldCME`(`AdjoinFieldClosure.lean`)が
   `Γ_{L_S} ≃ₜ* ↥S` を与え、Λ2 の `cyclotomeEquiv` がそれを円分子に運ぶ。
   最後に `S ≤ muFixer F (p^n)`(＝`μ_{p^n} ⊆ L_S`)から原始根を取り
   (`exists_isPrimitiveRoot_fixedField_of_le_muFixer`)、
   `μ_{p^n}(L_S)` が位数ちょうど `p^n` の巡回群であることを見て `ℤ/p^n` と同定する。

## ★★(ii) が残った理由(★測定結果。楽観的に書かない)

(ii) は「Artin 写像の**同変性**」である。`σ : L ≃ L` を `g ∈ Γ_F` の制限とすると、
必要なのは `Art_L(σ x) = Ψ_g(Art_L(x))`(`Ψ_g` は `g` による共役)であり、
Λ9 の `exists_abelianGalTorsion_equiv_cyclotomic` はこれを**仮説 `hΨ` として受け取っている**
(結論ではない)。★したがって Λ9 の到達点のうち、本ファイルが**そのまま使えたのは
仮説ゼロ版 `exists_abelianGalTorsion_equiv_rootsOfUnity` の方だけ**である。

★同変性が要る理由は「同型の取り方の不定性」ではない——それは Λ10 の
`pow_eq_pow_of_equivariant` が既に処理している(巡回群の自己同型群は可換なので、
`Λ_n(S)` の上で `g` が何乗として働くかは同型 `e` の取り方に依らない)。
残っているのは「**その冪指数が `χ_n(g)` である**」という、指数そのものの決定である。
これは Art の構成(素元 `π` と Lubin-Tate 級数 `f` の選択)を跨ぐ主張なので、
`σ(π)` が別の素元であること(`LubinTateUniformizerIndependence`)に加えて、
`σ` が形式群・捩れ点・不分岐部分を運ぶことが要る。★本ファイルの射程外である。

★**なお、この壁は原典の論拠(局所 Tate 双対性)とは独立である**——
原典の道なら `H^2(L, ℤ/p^n) ≅ ℤ/p^n` の不変写像が (ii) を直接与える。
すなわち (ii) は「局所相互律の同変性」か「局所双対性」のどちらかを通る。
2 つの道が同じ 1 点で止まっていることを本ファイルが確定させた。

## ★退化の自己検査

1. **`Λ_n(S)` は自明ではない**——`natCard_cyclotome_eq` が位数ちょうど `p^n` を与える。
   よって残った壁 (ii) は「自明群の上の作用」を主張する空虚な命題ではない。
2. **`S ≤ muFixer F (p^n)` は落とせない**。落とすと `μ_{p^n} ⊄ L_S` になり、
   位数は `p^n` を**割る**だけになる(Λ9 の `natCard_powTorsion_abelianGal_dvd`)。
   `K = ℚ_p`(`p` 奇)・`n = 1`・`S = ⊤` が反例(Λ10 のモジュール docstring)。
3. **`IsOpen S` は落とせない**。`fixedFieldLocalField` が開性を要求する
   (固定体が `ℚ_p` 上有限次であることに使う)。`S = ⊥` の退化は Λ10 の
   `not_nonempty_cyclotomeEquiv_bot` が排除している。
4. **`S.Normal` は (i) には要らない**——`nonempty_cyclotomeEquivZMod` は正規性を使わない。
   正規性が要るのは共役作用 `cyclotomeConj` を書く (ii) の側だけである。
   ★型に出ているとおりで、余分な仮定を混ぜていない。
5. **`n = 0` でも主張は成り立つ**(両辺とも自明群)。`p^0 = 1` なので `x^1 = 1 → x = 1`。
6. **壁の「易しい半分」は実際に証明してある**——
   `cyclotomeConj_eq_pow_cyclotomicCharacter_of_mem` が `g ∈ S` の場合を与える
   (両辺とも恒等)。★これは式の向き(`χ` か `χ^{-1}` か)を決めない——
   向きは Λ9 の `hΨ`(共変)に合わせてある。
7. **言い換えで強さを変えていない**——
   `cyclotomeConjIsCyclotomic_iff_torsionCyclotome` が Λ10 の壁との**同値性**を与える。

## 逸脱の記録

* 原典(pGC §1)は Proposition 1.1 の論拠を局所 Tate 双対性([2] Bloch-Kato Prop 3.8)に
  置いているが、本経路は**局所類体論**を通る。この逸脱は `ResearchPaper/pgc-goal.md` に
  記録済み(Λ10 と同じ)。
* `μ_{p^n}` の代わりに `Multiplicative (ZMod (p^n))` を使う点は Λ10 の逸脱記録 1 に従う。
  ★本ファイルはその同定を**実際に構成した**(Λ10 では仮定の中に隠れていた)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

open scoped IsMulCommutative

universe u v

/-! ## 1. 抽象核 —— 位相的アーベル化は「最大アーベル部分拡大の Galois 群」

★この節に局所体・付値・分岐の語彙は 1 つも出てこない。 -/

/-- **★抽象核**: `E/k` が(無限次でよい)Galois なら

  `Gal(E/k) ⧸ [Gal(E/k), Gal(E/k)]‾  ≅  Gal(E^{[G,G]‾}/k)`。

mathlib の `InfiniteGalois.normalAutEquivQuotient`(閉正規部分群による商と
固定体の Galois 群の一致)に `N := (commutator _).topologicalClosure` を代入したもの。
★**交換子群の位相的閉包**であること(単なる交換子群でないこと)がここで効く——
無限次 Galois 対応は**閉**部分群についてしか成り立たない。 -/
noncomputable def topAbelianizationEquivGalFixedField
    (k : Type u) (E : Type v) [Field k] [Field E] [Algebra k E] [IsGalois k E] :
    TopologicalAbelianization (E ≃ₐ[k] E) ≃*
      (IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure)
        ≃ₐ[k] IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure)) :=
  InfiniteGalois.normalAutEquivQuotient
    (⟨(commutator (E ≃ₐ[k] E)).topologicalClosure, Subgroup.isClosed_topologicalClosure _⟩)

/-- ★抽象核の値は**制限写像**である。

★次の節点((ii) の同変性)はここを使うことになる: 共役作用が `Gal(E/k)` の側で
どう見えるかは、制限写像の自然性で書けるからである。 -/
theorem topAbelianizationEquivGalFixedField_mk
    (k : Type u) (E : Type v) [Field k] [Field E] [Algebra k E] [IsGalois k E]
    (g : E ≃ₐ[k] E) :
    topAbelianizationEquivGalFixedField k E (QuotientGroup.mk g)
      = AlgEquiv.restrictNormalHom
          (IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure)) g :=
  InfiniteGalois.normalAutEquivQuotient_apply _ g

/-! ## 2. 具体層 —— `Γ_K^{ab}(位相的) ≅ Gal(K^ab/K)`

★Λ1/Λ2(群論側)と Λ8/Λ9(局所類体論側)を繋ぐ 1 本。 -/

variable {p : ℕ} [Fact p.Prime]

/-- **★★`Γ_K` の位相的アーベル化は `Gal(K^ab/K)`**。

`abelianClosure K` は定義から `fixedField ((commutator Γ_K).topologicalClosure)`
(`AbelianClosure.lean::abelianClosure_def`)なので、抽象核に `k := K`、`E := K̄` を
代入するだけで出る。★#59 を踏んでいない——中間体の層は 1 枚である。 -/
noncomputable def topAbelianizationEquivAbelianGal (K : PAdicLocalField p) :
    TopologicalAbelianization K.absGal ≃*
      (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) := by
  haveI := isGalois_closure K
  exact topAbelianizationEquivGalFixedField K.carrier K.closure

/-- 具体層でも値は制限写像。 -/
theorem topAbelianizationEquivAbelianGal_mk (K : PAdicLocalField p) (g : K.absGal) :
    topAbelianizationEquivAbelianGal K (QuotientGroup.mk g)
      = AlgEquiv.restrictNormalHom (abelianClosure K) g := by
  haveI := isGalois_closure K
  exact topAbelianizationEquivGalFixedField_mk K.carrier K.closure g

/-! ## 3. ★群論的円分子 `Λ_n(Γ_K)` は `μ_{p^n}(K)`

`cyclotome A m`(Λ2)と `powTorsion (A^{ab}) m`(Λ9)は**定義が同じ**
(どちらも `(powMonoidHom m).ker`)なので、橋は §2 の 1 本で足りる。 -/

/-- **★★★群論的に定義した円分子は `μ_{p^n}(K)` である**。

`cyclotome Γ_K (p^n)`(＝`tors_{p^n}(Γ_K^{ab})`、Λ2 が **`K` を一切見ずに**定義したもの)が
`μ_{p^n}(K)`(体の中の `p^n` 乗根)と同型であること。

★**これは同型の存在だけ**であって、`Γ`-作用の同変性は主張していない
(モジュール docstring の (ii))。★右辺は `K` に依存する群であり、`K = ℚ_p`(`p` 奇)なら
自明群である(Λ9 の退化の自己検査 1)。 -/
theorem nonempty_cyclotomeEquivRootsOfUnity (K : PAdicLocalField p) (n : ℕ) :
    Nonempty (↥(cyclotome K.absGal (p ^ n)) ≃* ↥(rootsOfUnity (p ^ n) K.carrier)) := by
  obtain ⟨e⟩ := exists_abelianGalTorsion_equiv_rootsOfUnity K n
  exact ⟨(powTorsionCongr (topAbelianizationEquivAbelianGal K) (p ^ n)).trans e⟩

/-- **開部分群版** —— `S ≤ Γ_K` が開なら `Λ_n(S) ≅ μ_{p^n}(L_S)`。

`absGalFixedFieldCME`(`AdjoinFieldClosure.lean`)が `Γ_{L_S} ≃ₜ* ↥S` を与える。
★`S.Normal` は要らない。 -/
theorem nonempty_cyclotomeSubgroupEquivRootsOfUnity (K : PAdicLocalField p)
    (S : Subgroup K.absGal) (hS : IsOpen (S : Set K.absGal)) (n : ℕ) :
    Nonempty (↥(cyclotome ↥S (p ^ n)) ≃*
      ↥(rootsOfUnity (p ^ n) (fixedFieldLocalField K S hS).carrier)) := by
  obtain ⟨e⟩ := nonempty_cyclotomeEquivRootsOfUnity (fixedFieldLocalField K S hS) n
  exact ⟨(cyclotomeEquiv (absGalFixedFieldCME K S hS) (p ^ n)).symm.trans e⟩

/-- **★★★(i) の到達点** —— `μ_{p^n} ⊆ L_S` なら `Λ_n(S) ≅ ℤ/p^n`。

`S ≤ muFixer K (p^n)` は「`L_S` が `p^n` 乗根をすべて含む」の言い換え
(Λ10 の `le_muFixer_iff`)なので、`L_S` の中に原始 `p^n` 乗根が取れて
`μ_{p^n}(L_S)` は位数ちょうど `p^n` の巡回群になる。

★これが Λ10 の `TorsionCyclotomeIsCyclotomic` の**前半**である。 -/
theorem nonempty_cyclotomeEquivZMod (K : PAdicLocalField p) {n : ℕ}
    {S : Subgroup K.absGal} (hS : IsOpen (S : Set K.absGal))
    (hle : S ≤ muFixer K (p ^ n)) :
    Nonempty (↥(cyclotome ↥S (p ^ n)) ≃* Multiplicative (ZMod (p ^ n))) := by
  haveI : NeZero (p ^ n) := ⟨pow_ne_zero _ (Fact.out : p.Prime).ne_zero⟩
  obtain ⟨η, hη⟩ := exists_isPrimitiveRoot_fixedField_of_le_muFixer K hS hle
  obtain ⟨e⟩ := nonempty_cyclotomeSubgroupEquivRootsOfUnity K S hS n
  haveI : IsCyclic ↥(rootsOfUnity (p ^ n) (fixedFieldLocalField K S hS).carrier) :=
    rootsOfUnity.isCyclic _ _
  refine ⟨e.trans (mulEquivOfCyclicCardEq ?_)⟩
  rw [Nat.card_eq_fintype_card, hη.card_rootsOfUnity, Nat.card_eq_fintype_card]
  simp

/-- ★**円分子の位数はちょうど `p^n`** —— 残った壁が空虚でないことの確認。 -/
theorem natCard_cyclotome_eq (K : PAdicLocalField p) {n : ℕ}
    {S : Subgroup K.absGal} (hS : IsOpen (S : Set K.absGal))
    (hle : S ≤ muFixer K (p ^ n)) :
    Nat.card ↥(cyclotome ↥S (p ^ n)) = p ^ n := by
  haveI : NeZero (p ^ n) := ⟨pow_ne_zero _ (Fact.out : p.Prime).ne_zero⟩
  obtain ⟨e⟩ := nonempty_cyclotomeEquivZMod K hS hle
  rw [Nat.card_congr e.toEquiv, Nat.card_eq_fintype_card]
  simp

/-! ## 4. ★★残った壁 —— 存在量化子を含まない 1 文 -/

def CyclotomeConjIsCyclotomic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★これだけが残った壁**(Λ10 の `TorsionCyclotomeIsCyclotomic` の後半)。

「`μ_{p^n} ⊆ L_S` なる開正規 `S ⊴ Γ_F` の円分子 `Λ_n(S) = tors_{p^n}(S^{ab})` の上で、
`g ∈ Γ_F` の共役作用は `x ↦ x^{χ_{F,n}(g)}` である。」

★Λ10 の `TorsionCyclotomeIsCyclotomic` と違い、**存在量化子を含まない**
(同型 `e` は `nonempty_cyclotomeEquivZMod` で構成済みなので、もう仮定に入れなくてよい)。

★古典的には Artin 写像の同変性 `Art_{L}(σ x) = σ Art_L(x) σ^{-1}` から出る。
Λ9 の `exists_abelianGalTorsion_equiv_cyclotomic` はこれを仮説 `hΨ` として受け取っており、
結論としては与えていない。★原典の道(局所 Tate 双対性)を通るなら、
`H^2(L, ℤ/p^n) ≅ ℤ/p^n` の不変写像がこれを直接与える。

★退化の自己検査: `S` の元 `g` に対しては両辺とも恒等
(左辺は `cyclotomeConj_coe_self`、右辺は `S ≤ muFixer F (p^n)` から `χ_n(g) = 1`)なので、
この主張は `Γ_F ⧸ S = Gal(L_S/F)` 上の主張である。 -/
def CyclotomeConjIsCyclotomic (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal),
    IsOpen (S : Set F.absGal) → S ≤ muFixer F (p ^ n) →
    ∀ (g : F.absGal) (x : ↥(cyclotome ↥S (p ^ n))),
      @cyclotomeConj F.absGal _ _ _ S hS (p ^ n) g x
        = x ^ ((PadicInt.toZModPow n
            ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val)

/-- **★Λ10 の仮定は (ii) だけに縮んだ**。

`TorsionCyclotomeIsCyclotomic`(同型 `e` の存在 + 作用が `χ` 倍)のうち、
同型の存在は `nonempty_cyclotomeEquivZMod` で構成済みなので、
残るのは作用の主張だけである。 -/
theorem torsionCyclotomeIsCyclotomic_of_cyclotomeConj
    (h : CyclotomeConjIsCyclotomic p) : TorsionCyclotomeIsCyclotomic p := by
  intro F n S hS hopen hle
  obtain ⟨e⟩ := nonempty_cyclotomeEquivZMod F hopen hle
  refine ⟨e, fun g x => ?_⟩
  rw [h F n S hS hopen hle g x, map_pow]

def cyclotomicCharacter_recoverable_of_cyclotomeConj.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★★★★★★★★★[pGC] Proposition 1.1 の現在地**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

円分指標の群論的復元は、いまや `CyclotomeConjIsCyclotomic p`(存在量化子を含まない 1 文)
だけに依存している。★`Skeleton/PGC/Section1.lean` の `sorry` は**まだ埋まっていない**——
埋めるにはこの 1 本が要る。 -/
theorem cyclotomicCharacter_recoverable_of_cyclotomeConj
    (h : CyclotomeConjIsCyclotomic p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicRecovery_of_torsionCyclotome (torsionCyclotomeIsCyclotomic_of_cyclotomeConj h)

/-- ★**逆向き**——Λ10 の `TorsionCyclotomeIsCyclotomic` からは (ii) が出る。

同型 `e` が作用を `c` 乗に写すなら、`e` は単射なので元の作用も `c` 乗である。 -/
theorem cyclotomeConjIsCyclotomic_of_torsionCyclotome
    (h : TorsionCyclotomeIsCyclotomic p) : CyclotomeConjIsCyclotomic p := by
  intro F n S hS hopen hle g x
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine e.injective ?_
  rw [he g x, map_pow]

/-- **★★★2 つの壁は同値**。

`CyclotomeConjIsCyclotomic`(存在量化子なし)と Λ10 の `TorsionCyclotomeIsCyclotomic`
(同型の存在を含む)は同値である。★すなわち本ファイルの言い換えで**主張は強くも弱くもなっていない**
——(i) の同型を実際に構成したぶんだけ、仮定から存在量化子が消えただけである。 -/
theorem cyclotomeConjIsCyclotomic_iff_torsionCyclotome :
    CyclotomeConjIsCyclotomic p ↔ TorsionCyclotomeIsCyclotomic p :=
  ⟨torsionCyclotomeIsCyclotomic_of_cyclotomeConj, cyclotomeConjIsCyclotomic_of_torsionCyclotome⟩

/-! ## 5. 退化の自己検査 -/

/-- ★**`muFixer` の上で `χ mod p^n` は `1`**。

`muFixer K (p^n)` は「`p^n` 乗根をすべて止める元」の集合だから、
原始根 `ζ` について `ζ = g ζ = ζ^{χ_n(g)}` となり `χ_n(g) = 1`。

★これは残った壁 `CyclotomeConjIsCyclotomic` の**無矛盾性検査**である:
`g ∈ S` のとき左辺(共役作用)は恒等(`cyclotomeConj_coe_self`)であり、
本補題により右辺の指数も `1` になる。★両辺が独立に恒等になることを確かめた。 -/
theorem toZModPow_cyclotomicCharacter_eq_one_of_mem_muFixer
    (K : PAdicLocalField p) (n : ℕ) {g : K.absGal} (hg : g ∈ muFixer K (p ^ n)) :
    (PadicInt.toZModPow n ((cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p]))) = 1 := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · haveI : Subsingleton (ZMod (p ^ 0)) := by rw [pow_zero]; infer_instance
    exact Subsingleton.elim _ _
  · set c : ZMod (p ^ n) :=
      PadicInt.toZModPow n ((cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p])) with hcdef
    obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_closure' K (p ^ n)
    have hgζ : g ζ = ζ := mem_muFixer.mp hg ζ hζ.pow_eq_one
    have hspec : g.toRingEquiv ζ = ζ ^ c.val :=
      cyclotomicCharacter.spec (n := n) p g.toRingEquiv ζ hζ.pow_eq_one
    have hpow : ζ ^ c.val = ζ ^ 1 := by
      rw [pow_one]
      exact hspec.symm.trans hgζ
    have hlt : (1 : ℕ) < p ^ n := Nat.one_lt_pow (by omega) (Fact.out : p.Prime).one_lt
    have hval : c.val = 1 := hζ.pow_inj (ZMod.val_lt c) hlt hpow
    have : ((c.val : ℕ) : ZMod (p ^ n)) = ((1 : ℕ) : ZMod (p ^ n)) := by rw [hval]
    rwa [ZMod.natCast_zmod_val, Nat.cast_one] at this

/-- ★★**残った壁は `S` の上では成り立っている**(壁の「易しい半分」)。

`g ∈ S` なら共役作用は恒等(`cyclotomeConj_coe_self`)、指数も `1`
(`toZModPow_cyclotomicCharacter_eq_one_of_mem_muFixer`)。★両辺が**独立の理由で**
一致するので、`CyclotomeConjIsCyclotomic` の式の向き(`χ` か `χ^{-1}` か)は
この検査では決まらない——★向きは Λ9 の `exists_abelianGalTorsion_equiv_cyclotomic` の
`hΨ`(共変、`Art(σ x) = σ Art(x) σ^{-1}`)に合わせてある。

★したがって残った壁の内容は `Γ_F ⧸ S = Gal(L_S/F)` の上にある。 -/
theorem cyclotomeConj_eq_pow_cyclotomicCharacter_of_mem
    (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal)
    (hle : S ≤ muFixer F (p ^ n)) {g : F.absGal} (hg : g ∈ S)
    (x : ↥(cyclotome ↥S (p ^ n))) :
    @cyclotomeConj F.absGal _ _ _ S hS (p ^ n) g x
      = x ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  have hself : @cyclotomeConj F.absGal _ _ _ S hS (p ^ n) g x = x :=
    cyclotomeConj_coe_self S (p ^ n) ⟨g, hg⟩ x
  rw [hself, toZModPow_cyclotomicCharacter_eq_one_of_mem_muFixer F n (hle hg)]
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · haveI : Subsingleton ↥(cyclotome (↥S) (p ^ 0)) := by
      constructor
      rintro ⟨a, ha⟩ ⟨b, hb⟩
      have ha' : a ^ (p ^ 0) = 1 := ha
      have hb' : b ^ (p ^ 0) = 1 := hb
      rw [pow_zero, pow_one] at ha' hb'
      exact Subtype.ext (ha'.trans hb'.symm)
    exact Subsingleton.elim _ _
  · haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
    have hlt : (1 : ℕ) < p ^ n := Nat.one_lt_pow (by omega) (Fact.out : p.Prime).one_lt
    rw [ZMod.val_one'' (by omega), pow_one]

end ABC3.Found.PGC

