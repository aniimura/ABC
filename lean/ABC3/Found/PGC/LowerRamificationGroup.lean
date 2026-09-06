import ABC3.Found.PGC.TotallyRamified
import Mathlib.RingTheory.Ideal.Over
import Mathlib.RingTheory.Ideal.Pointwise
import Mathlib.RingTheory.Nakayama
import Mathlib.RingTheory.DiscreteValuationRing.Basic

/-!
# 下付き分岐群 `G_n` と、完全分岐拡大の単項生成(Yoshida 2008 §5–§6)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108)。
本ファイルが形式化するのは

* **Lemma 5.11**(物理 p.12): `L'/L` が完全分岐で `α` が `L'` の素元なら `O_{L'} = O_L[α]`。
* **Definition 6.1**(物理 p.13–14): `i(σ) := v(σ(π) − π)`、`G_n := {σ | i(σ) > n}`。
  原典はさらに「`G_n` は正規で π の選び方に依らない。なぜなら
  `G_n = ker(G → Aut(O_{K'}/p^{n+1}))` だから」と述べる。

の 2 つ。

## ★逸脱の記録(1): `G_n` の定義を `Ideal.inertia` で与えた

原典は **π 依存**の定義 `G_n := {σ | v(σπ − π) > n}` から出発し、そのあとで
「正規性」と「π 非依存性」を、核としての特徴づけを経由して示す。

本ファイルは順序を入れ替え、**mathlib の `Ideal.inertia`**
(`Mathlib/RingTheory/Ideal/Defs.lean`、`I.inertia G = {σ | ∀ x, σ • x − x ∈ I}`)で

```
lowerRamificationGroup B G n := ((IsLocalRing.maximalIdeal B) ^ (n + 1)).inertia G
```

と定義し、原典の π による特徴づけを**定理側**(`mem_lowerRamificationGroup_iff` /
`mem_lowerRamificationGroup_iff_lt_addVal`)に回した。数学的な内容は同値である
(原典自身が `Def 6.1` の直後に上の核としての形を与えている)。

この入れ替えで消えるもの・残るものを実測した:

* **π 非依存性は完全に消える。** 定義に π が現れないので、二つの素元 `α`, `β` に
  よる特徴づけが一致することは `mem_pow_maximalIdeal_iff_of_adjoin_eq_top` として
  `.trans` 一発で出る(原典が `Prop 6.2` の証明で `σ(π')/π' = (σπ/π)(σu/u)` と
  計算している部分は、`G_n` の定義に関する限り不要になる)。
* **正規性は消えない。** ただし残るのは「環自己同型は極大イデアルを保つ」だけで
  (`smul_mem_maximalIdeal_pow`)、π の計算は要らない
  (`lowerRamificationGroup_normal`、4 行)。

## ★逸脱の記録(2): Lemma 5.11 を中村補題で証明した

原典の証明は「`v_{L'}(a_i α^i)` が相異なるので `v_{L'}(Σ a_i α^i) = min_i v_{L'}(a_i α^i)`」
という**付値の議論**で、`{1, α, …, α^{n−1}}` が基底であることを経由する。

本ファイルは同じ結論を**中村補題**で出した(`adjoin_uniformizer_eq_top`)。
`f = 1`(剰余体が伸びない)から `B = A[α] + 𝔪_B^k`(すべての `k`)を帰納法で得て
(`exists_mem_adjoin_sub_mem_pow`)、`𝔪_B^k ⊆ 𝔪_A B` となる `k` を取り
(`exists_pow_maximalIdeal_le_map`)、`Submodule.le_of_le_smul_of_le_jacobson_bot` を
当てる。`e = [L:K]` も基底の存在も使わないので、**仮定が原典より弱い**
(原典は `L'/L` が体の有限次拡大であることを使うが、こちらは
`Module.Finite A B` と `f = 1` と `𝔪_A B ≠ 0` だけ)。

★`ABC3.Found.Falt1.adjoin_eq_integralClosure_of_isEisensteinAt`
(`Found/Falt1/KaehlerAux.lean`)は**使わなかった**。あれは `PowerBasis K L` と
`minpoly R PB.gen` の Eisenstein 性を前提に取るので、「完全分岐 + 素元 ⇒ minpoly が
Eisenstein」を別途示す必要があり、遠回りになる(そちらの向きは
`Found/PGC/TotallyRamified.lean::isTotallyRamifiedAdjoin_of_eisenstein` が既にある)。

## ★危険信号への対応

* **完全分岐を落とすと `G_0 = ⊤` は偽**(`G_0` は惰性群になる)。本ファイルでは
  完全分岐は `Algebra.adjoin A {α} = ⊤`(= Lemma 5.11 の結論)として入っており、
  `lowerRamificationGroup_zero_eq_top` はそれを仮定に取る。
* **`Algebra.adjoin 𝒪_K {α} = 𝒪_L` を仮定として受け取らない。** これは
  `adjoin_uniformizer_eq_top` / `exists_uniformizer_adjoin_eq` の**結論**として供給する。
* **添字ずれ**: `G_n` は `𝔭^{n+1}` の inertia(`i(σ) > n ⟺ i(σ) ≥ n+1`)。
  `mem_pow_maximalIdeal_iff_lt_addVal` が `y ∈ 𝔪^{n+1} ↔ (n : ℕ∞) < addVal y` を
  与えて、原典の `i(σ) > n` と字面で一致させる。
* **退化検査**: `lowerRamificationGroup n := ⊤`(全 `n`)は
  `lowerRamificationGroup_zero_eq_top` と `lowerRamificationGroup_antitone` を通すが、
  `exists_lowerRamificationGroup_eq_bot` **だけ**を破る。
  ⇒ 6 番目を落とすと定義が退化する。

## 構成

抽象核(§1–§3)は「局所環 `B` と、`A`-代数として作用する群 `G`」だけで書いてある。
`PAdicLocalField` の重いインスタンス探索は §4 に閉じ込めてあり、抽象核の
`lean_check` は 0.05–0.6 秒で戻る(具体版は 0.5–1.8 秒)。
-/

namespace Ideal

/-- 惰性群はイデアルについて単調。 -/
theorem inertia_mono {B : Type*} [Ring B] {G : Type*} [Group G] [MulAction G B]
    {I J : Ideal B} (h : I ≤ J) : I.inertia G ≤ J.inertia G := by
  intro σ hσ
  rw [AddSubgroup.mem_inertia] at hσ ⊢
  exact fun x => h (hσ x)

end Ideal

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued

/-! ## §1 惰性群と単項生成 —— 抽象核

`Ideal.inertia` の membership は「すべての `x` について `σ • x − x ∈ I`」だが、
`B` が `A`-代数として `α` 一つで生成されていれば、`x = α` だけ見れば足りる。
これが原典 `Def 6.1` の「π による定義」と「核としての定義」を繋ぐ橋であり、
本ファイルで **Lemma 5.11 を使う唯一の場所**である。 -/

/-- `σ` が法 `I` で固定する元のなす `A`-部分代数。

`σ • (xy) − xy = (σ • x)(σ • y − y) + (σ • x − x) y` なので積で閉じ、
`σ` が `A`-線型なので `algebraMap A B` の像を含む。 -/
def fixedModSubalgebra (A : Type*) [CommRing A] {B : Type*} [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    (I : Ideal B) (σ : G) : Subalgebra A B where
  carrier := {x : B | σ • x - x ∈ I}
  mul_mem' {x y} hx hy := by
    simp only [Set.mem_setOf_eq, smul_mul'] at *
    have h : σ • x * σ • y - x * y = (σ • x) * (σ • y - y) + (σ • x - x) * y := by ring
    rw [h]
    exact I.add_mem (I.mul_mem_left _ hy) (I.mul_mem_right _ hx)
  add_mem' {x y} hx hy := by
    simp only [Set.mem_setOf_eq, smul_add] at *
    have h : σ • x + σ • y - (x + y) = (σ • x - x) + (σ • y - y) := by ring
    rw [h]; exact I.add_mem hx hy
  algebraMap_mem' r := by
    simp only [Set.mem_setOf_eq, Algebra.algebraMap_eq_smul_one, smul_comm σ r]
    simp

@[simp] theorem mem_fixedModSubalgebra {A : Type*} [CommRing A] {B : Type*} [CommRing B]
    [Algebra A B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {I : Ideal B} {σ : G} {x : B} :
    x ∈ fixedModSubalgebra A I σ ↔ σ • x - x ∈ I := Iff.rfl

/-- **`α` 一つで判定できる**——`B = A[α]` なら、惰性群の所属は `α` の像だけで決まる。
`←` は `fixedModSubalgebra` が `α` を含む部分代数であることから、
`→` は `α` を代入するだけ。 -/
theorem mem_inertia_iff_of_adjoin_eq_top {A : Type*} [CommRing A] {B : Type*} [CommRing B]
    [Algebra A B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    (I : Ideal B) {α : B} (hα : Algebra.adjoin A ({α} : Set B) = ⊤) (σ : G) :
    σ ∈ I.inertia G ↔ σ • α - α ∈ I := by
  constructor
  · intro h
    exact AddSubgroup.mem_inertia.mp h α
  · intro h
    rw [AddSubgroup.mem_inertia]
    intro x
    have hsub : (⊤ : Subalgebra A B) ≤ fixedModSubalgebra A I σ := by
      rw [← hα, Algebra.adjoin_le_iff]
      rintro y (rfl : y = α)
      exact h
    exact hsub Algebra.mem_top

/-! ## §2 Yoshida Lemma 5.11 —— 完全分岐拡大は素元 1 つで生成される

原典(物理 p.12)は
「If `L'/L` is totally ramified and `α` is a uniformizer of `L'`, then `O_{L'} = O_L[α]`」。
証明は本ファイル冒頭の「逸脱の記録(2)」の通り中村補題で書き換えてある。 -/

def exists_uniformizer_adjoin_eq.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 12, item := "Lemma 5.11", sectionId := "lemma-5-11" }

open IsLocalRing in
/-- **剰余体が伸びないときの逐次近似**——`x` は `A[α]` の元で `𝔪_B^k` を法として
近似できる(`k` について帰納法)。`𝔪_B^k = (α^k)` と `f = 1` を交互に使う。 -/
theorem exists_mem_adjoin_sub_mem_pow {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsLocalRing B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) :
    ∀ (k : ℕ) (x : B), ∃ c ∈ Algebra.adjoin A ({α} : Set B),
      x - c ∈ (maximalIdeal B) ^ k := by
  have hαmem : α ∈ maximalIdeal B := by rw [hα]; exact Ideal.mem_span_singleton_self α
  have hαadj : α ∈ Algebra.adjoin A ({α} : Set B) := Algebra.self_mem_adjoin_singleton A α
  intro k
  induction k with
  | zero => intro x; exact ⟨0, Subalgebra.zero_mem _, by simp⟩
  | succ k ih =>
    intro x
    obtain ⟨c, hc, hx⟩ := ih x
    rw [hα, Ideal.span_singleton_pow, Ideal.mem_span_singleton'] at hx
    obtain ⟨b, hb⟩ := hx
    obtain ⟨a, ha⟩ := hres b
    refine ⟨c + algebraMap A B a * α ^ k, ?_, ?_⟩
    · exact Subalgebra.add_mem _ hc
        (Subalgebra.mul_mem _ (Subalgebra.algebraMap_mem _ a) (Subalgebra.pow_mem _ hαadj k))
    · have key : x - (c + algebraMap A B a * α ^ k) = α ^ k * (b - algebraMap A B a) := by
        have hx' : x = b * α ^ k + c := by rw [hb]; ring
        rw [hx']; ring
      rw [key, pow_succ]
      exact Ideal.mul_mem_mul (Ideal.pow_mem_pow hαmem k) ha

open IsLocalRing in
/-- **Yoshida Lemma 5.11(抽象核)**——`B = A[α]`。

仮定は「`A` 局所・`B` 局所・`B` は `A` 上有限・`α` は `𝔪_B` の生成元・
剰余体が伸びない(`hres`)・`𝔪_A B` のある冪より上が取れる(`hpow`)」だけ。
`Submodule.le_of_le_smul_of_le_jacobson_bot`(中村)に
`N := A[α]`、`N' := ⊤`、`I := 𝔪_A` を当てる。 -/
theorem adjoin_uniformizer_eq_top {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B]
    [Algebra A B] [IsLocalRing B] [Module.Finite A B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α})
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hpow : ∃ k : ℕ, (maximalIdeal B) ^ k ≤ Ideal.map (algebraMap A B) (maximalIdeal A)) :
    Algebra.adjoin A ({α} : Set B) = ⊤ := by
  obtain ⟨k, hk⟩ := hpow
  have hjac : (maximalIdeal A) ≤ (⊥ : Ideal A).jacobson := by
    rw [IsLocalRing.jacobson_eq_maximalIdeal (⊥ : Ideal A) bot_ne_top]
  have hNN : (⊤ : Submodule A B) ≤
      Subalgebra.toSubmodule (Algebra.adjoin A ({α} : Set B))
        ⊔ (maximalIdeal A) • (⊤ : Submodule A B) := by
    intro x _
    obtain ⟨c, hc, hx⟩ := exists_mem_adjoin_sub_mem_pow hα hres k x
    have hxc : x = c + (x - c) := by ring
    rw [hxc]
    refine Submodule.add_mem_sup ((Subalgebra.mem_toSubmodule _).mpr hc) ?_
    rw [Ideal.smul_top_eq_map]
    exact hk hx
  have hle := Submodule.le_of_le_smul_of_le_jacobson_bot Module.Finite.fg_top hjac hNN
  exact eq_top_iff.mpr fun x _ => (Subalgebra.mem_toSubmodule _).mp (hle Submodule.mem_top)

open IsLocalRing in
/-- DVR では、`𝔪_A` の像が `0` でなければ `𝔪_B` のある冪がその中に入る
(`adjoin_uniformizer_eq_top` の `hpow` を供給する)。 -/
theorem exists_pow_maximalIdeal_le_map {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B]
    [IsDomain B] [IsDiscreteValuationRing B] [Algebra A B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α})
    (hne : Ideal.map (algebraMap A B) (maximalIdeal A) ≠ ⊥) :
    ∃ k : ℕ, (maximalIdeal B) ^ k ≤ Ideal.map (algebraMap A B) (maximalIdeal A) := by
  obtain ⟨n, hn⟩ := IsDiscreteValuationRing.ideal_eq_span_pow_irreducible hne
    ((IsDiscreteValuationRing.irreducible_iff_uniformizer α).mpr hα)
  exact ⟨n, by rw [hα, Ideal.span_singleton_pow, hn]⟩

open IsLocalRing in
/-- **Yoshida Lemma 5.11(存在形)**——素元 `α` が取れて `B = A[α]`。

★`Algebra.adjoin A {α} = ⊤` を**仮定ではなく結論として**供給するのがこの補題の役目。 -/
theorem exists_uniformizer_adjoin_eq {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B]
    [IsDomain B] [IsDiscreteValuationRing B] [Algebra A B] [Module.Finite A B]
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hne : Ideal.map (algebraMap A B) (maximalIdeal A) ≠ ⊥) :
    ∃ α : B, maximalIdeal B = Ideal.span {α} ∧ Algebra.adjoin A ({α} : Set B) = ⊤ := by
  obtain ⟨α, hα⟩ := IsDiscreteValuationRing.exists_irreducible B
  have hspan : maximalIdeal B = Ideal.span {α} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer α).mp hα
  exact ⟨α, hspan, adjoin_uniformizer_eq_top hspan hres (exists_pow_maximalIdeal_le_map hspan hne)⟩

/-! ## §3 Yoshida Definition 6.1 —— 下付き分岐群 `G_n`

原典(物理 p.13–14)は `i(σ) := v(σ(π) − π)`(`i(id) = ∞`)と置き、`n ≥ 0` に対して
`G_n := {σ ∈ G | i(σ) > n}` と定義する。続けて、これらが正規で π の選び方に依らない
ことを「`G_n = {σ | v(σ(a) − a) > n for all a ∈ O_{K'}}` が
`G → Aut(O_{K'}/p^{n+1}_{K'})` の核だから」として示す。 -/

def lowerRamificationGroup.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Definition 6.1", sectionId := "def-6-1" }

open IsLocalRing in
/-- **下付き分岐群** `G_n := ((𝔪_B)^{n+1}).inertia G`。

★★逸脱: 原典は π 依存の定義から出発する。ここでは原典自身が与えている
「核としての形」を定義に採り、π による特徴づけを
`mem_lowerRamificationGroup_iff` に回した(ファイル冒頭「逸脱の記録(1)」)。
★添字: `G_n` は `𝔪^{n+1}` の inertia(`i(σ) > n ⟺ i(σ) ≥ n+1`)。 -/
def lowerRamificationGroup (B : Type*) [CommRing B] [IsLocalRing B] (G : Type*) [Group G]
    [MulSemiringAction G B] (n : ℕ) : Subgroup G :=
  ((maximalIdeal B) ^ (n + 1)).inertia G

open IsLocalRing in
@[simp] theorem mem_lowerRamificationGroup_iff_forall {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {n : ℕ} {σ : G} :
    σ ∈ lowerRamificationGroup B G n ↔ ∀ x : B, σ • x - x ∈ (maximalIdeal B) ^ (n + 1) :=
  AddSubgroup.mem_inertia

open IsLocalRing in
/-- **Yoshida Definition 6.1 の橋**——`B = A[α]` のとき `σ ∈ G_n ⟺ σα − α ∈ 𝔪^{n+1}`。
★ここでだけ §2(Lemma 5.11)の結論を使う。 -/
theorem mem_lowerRamificationGroup_iff {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : Algebra.adjoin A ({α} : Set B) = ⊤) (n : ℕ) (σ : G) :
    σ ∈ lowerRamificationGroup B G n ↔ σ • α - α ∈ (maximalIdeal B) ^ (n + 1) :=
  mem_inertia_iff_of_adjoin_eq_top _ hα σ

open IsLocalRing IsDiscreteValuationRing in
/-- `y ∈ 𝔪^{n+1} ⟺ n < v(y)`——原典の `i(σ) > n` の字面に合わせるための翻訳。 -/
theorem mem_pow_maximalIdeal_iff_lt_addVal {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (n : ℕ) (y : B) :
    y ∈ (maximalIdeal B) ^ (n + 1) ↔ (n : ℕ∞) < addVal B y := by
  have hirr : Irreducible α := (IsDiscreteValuationRing.irreducible_iff_uniformizer α).mpr hα
  rw [hα, Ideal.span_singleton_pow, Ideal.mem_span_singleton, ← addVal_le_iff_dvd,
    addVal_pow, addVal_uniformizer hirr, nsmul_eq_mul, mul_one, Nat.cast_add, Nat.cast_one,
    ENat.add_one_le_iff (by simp)]

open IsLocalRing IsDiscreteValuationRing in
/-- **Yoshida Definition 6.1 そのままの形**: `σ ∈ G_n ⟺ i(σ) > n`
(`i(σ) = v(σα − α)`、`i(id) = ∞` は `addVal B 0 = ⊤` が担う)。 -/
theorem mem_lowerRamificationGroup_iff_lt_addVal {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hα : Algebra.adjoin A ({α} : Set B) = ⊤)
    (n : ℕ) (σ : G) :
    σ ∈ lowerRamificationGroup B G n ↔ (n : ℕ∞) < addVal B (σ • α - α) :=
  (mem_lowerRamificationGroup_iff hα n σ).trans
    (mem_pow_maximalIdeal_iff_lt_addVal huni n _)

open IsLocalRing in
/-- **π 非依存性**——原典 `Def 6.1` が核としての形を経由して示している部分。
`Ideal.inertia` による定義には π が現れないので `.trans` 一発で出る。 -/
theorem mem_pow_maximalIdeal_iff_of_adjoin_eq_top {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α β : B} (hα : Algebra.adjoin A ({α} : Set B) = ⊤)
    (hβ : Algebra.adjoin A ({β} : Set B) = ⊤) (n : ℕ) (σ : G) :
    σ • α - α ∈ (maximalIdeal B) ^ (n + 1) ↔ σ • β - β ∈ (maximalIdeal B) ^ (n + 1) :=
  (mem_lowerRamificationGroup_iff hα n σ).symm.trans (mem_lowerRamificationGroup_iff hβ n σ)

open IsLocalRing in
/-- 局所環の環自己同型は極大イデアルを保つ(単元を単元に写すから)。 -/
theorem smul_mem_maximalIdeal {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] (σ : G) {x : B} (hx : x ∈ maximalIdeal B) :
    σ • x ∈ maximalIdeal B := by
  rw [IsLocalRing.mem_maximalIdeal] at hx ⊢
  intro hu
  exact hx (by simpa using hu.map (MulSemiringAction.toRingHom G B σ⁻¹))

open IsLocalRing in
/-- 極大イデアルの冪も保つ。 -/
theorem smul_mem_maximalIdeal_pow {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] (σ : G) (n : ℕ) {x : B} (hx : x ∈ (maximalIdeal B) ^ n) :
    σ • x ∈ (maximalIdeal B) ^ n := by
  induction n generalizing x with
  | zero => simp
  | succ n ih =>
    rw [pow_succ] at hx ⊢
    refine Submodule.smul_induction_on hx ?_ ?_
    · intro a ha b hb
      rw [smul_eq_mul, smul_mul']
      exact Ideal.mul_mem_mul (ih ha) (smul_mem_maximalIdeal σ hb)
    · intro y z hy hz; rw [smul_add]; exact Ideal.add_mem _ hy hz

/-- **`G_n` は `G` の正規部分群**。★原典は π 依存の定義から出発するので
`σ(π')/π' = (σπ/π)(σu/u)` の計算が要るが、`Ideal.inertia` で定義すると
「自己同型が `𝔪` を保つ」だけで済む(π は現れない)。 -/
instance lowerRamificationGroup_normal (B : Type*) [CommRing B] [IsLocalRing B] (G : Type*)
    [Group G] [MulSemiringAction G B] (n : ℕ) : (lowerRamificationGroup B G n).Normal where
  conj_mem σ hσ τ := by
    rw [mem_lowerRamificationGroup_iff_forall] at hσ ⊢
    intro y
    have h1 : (τ * σ * τ⁻¹) • y - y = τ • (σ • (τ⁻¹ • y) - τ⁻¹ • y) := by
      simp [mul_smul, smul_sub]
    rw [h1]
    exact smul_mem_maximalIdeal_pow τ (n + 1) (hσ (τ⁻¹ • y))

open IsLocalRing in
/-- **`G = G_0`(完全分岐のとき)**。原典 `Def 6.1` の
「Then `G = G_0` as `K'/K` is totally ramified」。

★完全分岐性は `hα : Algebra.adjoin A {α} = ⊤`(= Lemma 5.11 の結論)として
入っている。これを落とすと `G_0` は惰性群になり主張は**偽**になる。 -/
theorem lowerRamificationGroup_zero_eq_top {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hαm : α ∈ maximalIdeal B) (hα : Algebra.adjoin A ({α} : Set B) = ⊤) :
    lowerRamificationGroup B G 0 = ⊤ := by
  rw [eq_top_iff]
  intro σ _
  rw [mem_lowerRamificationGroup_iff hα]
  simpa using Ideal.sub_mem _ (smul_mem_maximalIdeal σ hαm) hαm

/-- `α` が `B` を生成し作用が忠実なら、`σ • α = α` から `σ = 1`。 -/
theorem eq_one_of_smul_eq_of_adjoin_eq_top {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B]
    {α : B} (hα : Algebra.adjoin A ({α} : Set B) = ⊤) {σ : G} (h : σ • α = α) : σ = 1 := by
  have hmem : σ ∈ (⊥ : Ideal B).inertia G :=
    (mem_inertia_iff_of_adjoin_eq_top ⊥ hα σ).mpr (by rw [h]; simp)
  rw [AddSubgroup.mem_inertia] at hmem
  refine eq_of_smul_eq_smul (α := B) fun x => ?_
  have hx : σ • x - x ∈ (⊥ : Ideal B) := hmem x
  rw [Ideal.mem_bot, sub_eq_zero] at hx
  simpa using hx

/-- **`G_n` は減少列**(`Antitone`)。 -/
theorem lowerRamificationGroup_antitone (B : Type*) [CommRing B] [IsLocalRing B] (G : Type*)
    [Group G] [MulSemiringAction G B] : Antitone (lowerRamificationGroup B G) := by
  intro m n hmn
  exact Ideal.inertia_mono (Ideal.pow_le_pow_right (Nat.succ_le_succ hmn))

open IsLocalRing in
/-- **`G_n = {id}` for sufficiently large `n`**(原典 `Def 6.1`)。

★退化検査: `lowerRamificationGroup n := ⊤`(全 `n`)は
`lowerRamificationGroup_zero_eq_top` と `lowerRamificationGroup_antitone` を通すが
**この定理だけ**を破る。⇒ これを落とすと定義が退化する。

証明は Krull の交叉定理 `⨅ i, 𝔪^i = ⊥`(`B` は Noether 局所環)と、
`G` が有限であることから。 -/
theorem exists_lowerRamificationGroup_eq_bot {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsLocalRing B] [IsNoetherianRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    {α : B} (hα : Algebra.adjoin A ({α} : Set B) = ⊤) :
    ∃ N : ℕ, ∀ n, N ≤ n → lowerRamificationGroup B G n = ⊥ := by
  classical
  have hiInf : (⨅ i : ℕ, (maximalIdeal B) ^ i) = ⊥ :=
    Ideal.iInf_pow_eq_bot_of_isLocalRing _ (IsLocalRing.maximalIdeal.isMaximal B).ne_top
  have key : ∀ σ : G, ∃ n : ℕ, ∀ m, n ≤ m → σ ∈ lowerRamificationGroup B G m → σ = 1 := by
    intro σ
    by_cases h1 : σ = 1
    · exact ⟨0, fun _ _ _ => h1⟩
    · have hnz : σ • α - α ≠ 0 := fun h =>
        h1 (eq_one_of_smul_eq_of_adjoin_eq_top hα (by rwa [sub_eq_zero] at h))
      have hex : ∃ i : ℕ, σ • α - α ∉ (maximalIdeal B) ^ i := by
        by_contra hc
        rw [not_exists] at hc
        have hmem : σ • α - α ∈ ⨅ i : ℕ, (maximalIdeal B) ^ i :=
          Ideal.mem_iInf.mpr fun i => not_not.mp (hc i)
        rw [hiInf, Ideal.mem_bot] at hmem
        exact hnz hmem
      obtain ⟨i, hi⟩ := hex
      refine ⟨i, fun m hm hmem => absurd ?_ hi⟩
      exact Ideal.pow_le_pow_right (le_trans hm (Nat.le_succ m))
        ((mem_lowerRamificationGroup_iff hα m σ).mp hmem)
  choose f hf using key
  haveI : Fintype G := Fintype.ofFinite G
  refine ⟨Finset.sup Finset.univ f, fun n hn => ?_⟩
  rw [eq_bot_iff]
  intro σ hσ
  exact Subgroup.mem_bot.mpr
    (hf σ n (le_trans (Finset.le_sup (Finset.mem_univ σ)) hn) hσ)

open IsLocalRing in
/-- **部分群への制限**——`H ≤ G` に対し `(G_n).subgroupOf H = H_n`。★`rfl`。
Sen の `H_n := G_n ∩ ⟨σ⟩`(Hasse–Arf の心臓)がこの形を要求する。 -/
theorem subgroupOf_lowerRamificationGroup {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] (H : Subgroup G) (n : ℕ) :
    (lowerRamificationGroup B G n).subgroupOf H = lowerRamificationGroup B H n := rfl

/-! ## §4 `PAdicLocalField` への具体化

`A := 𝒪[K.carrier]`、`B := adjoinIntegers K x`(`= 𝒪_{K(x)}`)、
`G := Gal(K(x)/K)` に当てはめる。既存在庫:

* `isDiscreteValuationRing_carrierIntegers` / `isDiscreteValuationRing_adjoinIntegers`
* `instIsLocalRingAdjoinIntegers` / `adjoinIntegersAlgebra` / `module_finite_adjoinIntegers`
* `isIntegral_iff_norm_le_one`(`Found/PGC/UnramifiedExtension.lean`)

新しく作るのは **`Gal(K(x)/K)` の `adjoinIntegers K x` への環作用**
(`mulSemiringActionAdjoinIntegers`)と、完全分岐からの剰余体全射性。 -/

variable {p : ℕ} [Fact p.Prime]

/-- **完全分岐なら剰余体は伸びない**——`𝒪_K → 𝒪_L/𝔭_L` が全射。
`inertiaDegree K x = 1` から剰余体の濃度が等しく、単射な体の射は全単射。 -/
theorem exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) (b : adjoinIntegers K x) :
    ∃ a : 𝒪[K.carrier], b - algebraMap 𝒪[K.carrier] (adjoinIntegers K x) a
      ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) := by
  classical
  haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K x)) := Fintype.ofFinite _
  set f : 𝓀[K.carrier] →+* IsLocalRing.ResidueField (adjoinIntegers K x) :=
    IsLocalRing.ResidueField.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) with hf
  have hcard : Fintype.card 𝓀[K.carrier]
      = Fintype.card (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
    have h1 := (isTotallyRamifiedAdjoin_iff_residueDegree K x).mp ht
    rw [residueDegree] at h1
    rw [← Nat.card_eq_fintype_card, ← Nat.card_eq_fintype_card, h1]
  have hbij : Function.Bijective f :=
    (Fintype.bijective_iff_injective_and_card f).mpr ⟨f.injective, hcard⟩
  obtain ⟨c, hc⟩ := hbij.2 (IsLocalRing.residue _ b)
  obtain ⟨a, rfl⟩ := IsLocalRing.residue_surjective (R := 𝒪[K.carrier]) c
  refine ⟨a, ?_⟩
  rw [← IsLocalRing.residue_eq_zero_iff, map_sub, sub_eq_zero]
  rw [hf, IsLocalRing.ResidueField.map_residue] at hc
  exact hc.symm

/-- `𝒪[K.carrier] → adjoinIntegers K x` は単射。 -/
theorem injective_algebraMap_adjoinIntegers (K : PAdicLocalField p) (x : K.closure) :
    Function.Injective (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) := by
  intro a b hab
  have h : algebraMap K.carrier K.closure (a : K.carrier)
      = algebraMap K.carrier K.closure (b : K.carrier) :=
    congrArg (fun z => ((z : adjoinIntegers K x) : K.closure)) hab
  exact Subtype.ext ((algebraMap K.carrier K.closure).injective h)

/-- `𝔪_K` の像は `0` でない(`adjoin_uniformizer_eq_top` の `hpow` を出すのに要る)。 -/
theorem map_maximalIdeal_ne_bot_adjoinIntegers (K : PAdicLocalField p) (x : K.closure) :
    Ideal.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (IsLocalRing.maximalIdeal 𝒪[K.carrier]) ≠ ⊥ := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  obtain ⟨π, hπ⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  have hπm : π ∈ IsLocalRing.maximalIdeal 𝒪[K.carrier] :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer π).mp hπ ▸
      Ideal.mem_span_singleton_self π
  intro hc
  have hmem : algebraMap 𝒪[K.carrier] (adjoinIntegers K x) π
      ∈ Ideal.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
        (IsLocalRing.maximalIdeal 𝒪[K.carrier]) := Ideal.mem_map_of_mem _ hπm
  rw [hc, Ideal.mem_bot] at hmem
  exact hπ.ne_zero (injective_algebraMap_adjoinIntegers K x (by simpa using hmem))

/-- **Yoshida Lemma 5.11 の `PAdicLocalField` 版**——
`K(x)/K` 完全分岐、`α` が `𝒪_{K(x)}` の素元 ⇒ `𝒪_K[α] = 𝒪_{K(x)}`。 -/
theorem adjoin_uniformizer_eq_top_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {α : adjoinIntegers K x}
    (hα : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α}) :
    Algebra.adjoin 𝒪[K.carrier] ({α} : Set (adjoinIntegers K x)) = ⊤ := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI := module_finite_adjoinIntegers K x
  haveI := isDiscreteValuationRing_carrierIntegers K
  exact adjoin_uniformizer_eq_top hα
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht)
    (exists_pow_maximalIdeal_le_map hα (map_maximalIdeal_ne_bot_adjoinIntegers K x))

/-! ### `Gal(K(x)/K)` の `𝒪_{K(x)}` への作用 -/

theorem mem_adjoinIntegers_iff' (K : PAdicLocalField p) (x : K.closure)
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    y ∈ adjoinIntegers K x ↔ ‖y‖ ≤ 1 := Iff.rfl

/-- `K(x)` の `K`-自己同型は整数環を保つ——`‖y‖ ≤ 1 ⟺ y が 𝒪_K 上整`
(`isIntegral_iff_norm_le_one`)を経由すると、σ が `𝒪_K`-代数射であることだけで済む。 -/
theorem smul_mem_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    {y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)}
    (hy : y ∈ adjoinIntegers K x) : σ y ∈ adjoinIntegers K x := by
  haveI : IsScalarTower 𝒪[K.carrier] K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  rw [mem_adjoinIntegers_iff', ← isIntegral_iff_norm_le_one] at hy ⊢
  exact hy.map (σ.toAlgHom.restrictScalars 𝒪[K.carrier])

/-- **`Gal(K(x)/K)` は `𝒪_{K(x)}` に環として作用する。** -/
noncomputable instance mulSemiringActionAdjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    MulSemiringAction
      ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (adjoinIntegers K x) where
  smul σ y := ⟨σ (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
    smul_mem_adjoinIntegers K x σ y.2⟩
  one_smul _ := Subtype.ext rfl
  mul_smul _ _ _ := Subtype.ext rfl
  smul_zero σ := Subtype.ext (map_zero σ)
  smul_add σ _ _ := Subtype.ext (map_add σ _ _)
  smul_one σ := Subtype.ext (map_one σ)
  smul_mul σ _ _ := Subtype.ext (map_mul σ _ _)

@[simp] theorem coe_smul_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (y : adjoinIntegers K x) :
    ((σ • y : adjoinIntegers K x) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = σ (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := rfl

/-- `Gal(K(x)/K)` の作用は `𝒪_K`-線型。 -/
instance smulCommClassAdjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    SMulCommClass
      ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      𝒪[K.carrier] (adjoinIntegers K x) where
  smul_comm σ a y := by
    apply Subtype.ext
    rw [Algebra.smul_def, Algebra.smul_def]
    show σ (((algebraMap 𝒪[K.carrier] (adjoinIntegers K x) a : adjoinIntegers K x) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        * (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      = ((algebraMap 𝒪[K.carrier] (adjoinIntegers K x) a : adjoinIntegers K x) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) * σ (y : _)
    rw [map_mul]
    exact congrArg
      (fun z : IntermediateField.adjoin K.carrier ({x} : Set K.closure) =>
        z * σ (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (σ.commutes (a : K.carrier))

/-- `Gal(K(x)/K)` の `𝒪_{K(x)}` への作用は忠実——`𝒪_{K(x)}` は `K(x)` の付値環なので
`y` か `y⁻¹` の一方は必ず入る。 -/
instance faithfulSMulAdjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    FaithfulSMul
      ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (adjoinIntegers K x) where
  eq_of_smul_eq_smul {σ τ} h := by
    have hint : ∀ y : IntermediateField.adjoin K.carrier ({x} : Set K.closure),
        ‖y‖ ≤ 1 → σ y = τ y := by
      intro y hy
      exact congrArg
        (fun z : adjoinIntegers K x =>
          (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
        (h ⟨y, hy⟩)
    refine AlgEquiv.ext fun y => ?_
    by_cases hy : ‖y‖ ≤ 1
    · exact hint y hy
    · have hinv : ‖y⁻¹‖ ≤ 1 := by
        rw [norm_inv]
        exact inv_le_one_of_one_le₀ (le_of_lt (not_le.mp hy))
      have hh := hint y⁻¹ hinv
      rw [map_inv₀, map_inv₀, inv_inj] at hh
      exact hh

/-! ### 下付き分岐群(`PAdicLocalField` 版) -/

/-- **`PAdicLocalField` 版の下付き分岐群** `G_n ⊆ Gal(K(x)/K)`。 -/
noncomputable def lowerRamificationGroupAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (n : ℕ) : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  lowerRamificationGroup (adjoinIntegers K x) _ n

/-- **Yoshida Definition 6.1(`PAdicLocalField` 版)**。 -/
theorem mem_lowerRamificationGroupAdjoin_iff (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {α : adjoinIntegers K x}
    (hα : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α}) (n : ℕ)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    σ ∈ lowerRamificationGroupAdjoin K x n ↔
      σ • α - α ∈ (IsLocalRing.maximalIdeal (adjoinIntegers K x)) ^ (n + 1) :=
  mem_lowerRamificationGroup_iff (adjoin_uniformizer_eq_top_adjoinIntegers K x ht hα) n σ

/-- **`G_0 = ⊤`(完全分岐)**。 -/
theorem lowerRamificationGroupAdjoin_zero_eq_top (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) : lowerRamificationGroupAdjoin K x 0 = ⊤ := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  obtain ⟨α, hα⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K x)
  have hspan : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer α).mp hα
  exact lowerRamificationGroup_zero_eq_top (hspan ▸ Ideal.mem_span_singleton_self α)
    (adjoin_uniformizer_eq_top_adjoinIntegers K x ht hspan)

/-- **`G_n` は減少列**。 -/
theorem lowerRamificationGroupAdjoin_antitone (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Antitone (lowerRamificationGroupAdjoin K x) :=
  lowerRamificationGroup_antitone _ _

/-- **`G_n` は十分大きな `n` で自明**。 -/
theorem exists_lowerRamificationGroupAdjoin_eq_bot (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ∃ N : ℕ, ∀ n, N ≤ n → lowerRamificationGroupAdjoin K x n = ⊥ := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI : IsNoetherianRing (adjoinIntegers K x) := inferInstance
  obtain ⟨α, hα⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K x)
  have hspan : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer α).mp hα
  exact exists_lowerRamificationGroup_eq_bot
    (adjoin_uniformizer_eq_top_adjoinIntegers K x ht hspan)

/-- **部分群への制限**(Sen の `H_n := G_n ∩ ⟨σ⟩` が要る形)。★`rfl`。 -/
theorem subgroupOf_lowerRamificationGroupAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (H : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) (n : ℕ) :
    (lowerRamificationGroupAdjoin K x n).subgroupOf H
      = lowerRamificationGroup (adjoinIntegers K x) H n := rfl

/-! ### 付値の形(原典の字面)

`IsDiscreteValuationRing.addVal` を statement に出すには DVR インスタンスが要るので、
ここだけ `attribute [local instance]` で入れる。 -/

attribute [local instance] isDiscreteValuationRing_adjoinIntegers
attribute [local instance] isDiscreteValuationRing_carrierIntegers
attribute [local instance] module_finite_adjoinIntegers

/-- **Yoshida Lemma 5.11 の `PAdicLocalField` 版(存在形)**——
`K(x)/K` が完全分岐なら、素元 `α` が取れて `𝒪_K[α] = 𝒪_{K(x)}`。 -/
theorem exists_uniformizer_adjoin_eq_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ∃ α : adjoinIntegers K x,
      IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α} ∧
      Algebra.adjoin 𝒪[K.carrier] ({α} : Set (adjoinIntegers K x)) = ⊤ :=
  exists_uniformizer_adjoin_eq
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht)
    (map_maximalIdeal_ne_bot_adjoinIntegers K x)

/-- **Yoshida Definition 6.1 の `PAdicLocalField` 版(付値の形)**:
`σ ∈ G_n ⟺ n < v_{K(x)}(σα − α)`。 -/
theorem mem_lowerRamificationGroupAdjoin_iff_lt_addVal (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {α : adjoinIntegers K x}
    (hα : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α}) (n : ℕ)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    σ ∈ lowerRamificationGroupAdjoin K x n ↔
      (n : ℕ∞) < IsDiscreteValuationRing.addVal (adjoinIntegers K x) (σ • α - α) :=
  mem_lowerRamificationGroup_iff_lt_addVal hα
    (adjoin_uniformizer_eq_top_adjoinIntegers K x ht hα) n σ

end ABC3.Found.PGC
