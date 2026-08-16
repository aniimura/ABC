import ABC3.Found.FrdI.MonoidPrime
import ABC3.Found.FrdI.Frobenioid

/-!
# [FrdI] Definition 2.4 —— perf-factorial / `Λ supports M` / `C(d)`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.47–p.48。

原文 (FrdI p.47):
> Definition 2.4.

## ★段 1 —— 実数化 `⊗ R` の実体

★★**原文は `M rlf` を `M pf ⊗R` の一言で済ませている**(p.47 の (c))。
★**その一言の中身をここで作る。**

★`M^pf` の元は形式的な分数 `m/n`(`Pf.mk m n`)なので、
`f : M →+ ℝ≥0` が与えられれば `M^pf → ℝ≥0` は **`m/n ↦ f m / n`** で延びる
(`pfLift`)。★**これが `⊗ ℝ≥0` の計算的な実体**である。

★★**`M` が monoprime なら `M ≃+ Λ≥0` であり、`Λ≥0` はどの場合も
`ℝ≥0` に自然に埋まる**ので、`M^pf ↪ ℝ≥0` が得られる(`monoprimeLift`)。

★★★**同型の取り方に依らないこと**は `pfLift_unique_of_smul` で押さえる ——
2 つの埋め込みは **`ℝ≥0` の正の実数倍だけ違う**。
-/

namespace ABC3.Found.FrdI

open scoped NNReal
open CategoryTheory

universe v u w u2 v2

variable (M : Type w) [AddCommMonoid M]

/-! ## ★`M^pf → ℝ≥0` への延長 -/

variable {M}

/-- ★`M →+ ℝ≥0` を `M^pf → ℝ≥0` へ延ばす(台の写像)。 -/
noncomputable def pfLiftFun (f : M →+ ℝ≥0) : Pf M → ℝ≥0 :=
  Quotient.lift (fun x : M × ℕ+ => f x.1 / ((x.2 : ℕ) : ℝ≥0)) (by
    rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨k, hk⟩
    have hk0 : ((k : ℕ) : ℝ≥0) ≠ 0 := by
      simpa using (k.pos.ne' : (0 : ℕ) ≠ (k : ℕ))
    have ha0 : ((a : ℕ) : ℝ≥0) ≠ 0 := by
      simpa using (a.pos.ne' : (0 : ℕ) ≠ (a : ℕ))
    have ha'0 : ((a' : ℕ) : ℝ≥0) ≠ 0 := by
      simpa using (a'.pos.ne' : (0 : ℕ) ≠ (a' : ℕ))
    have hf := congrArg f hk
    rw [map_nsmul, map_nsmul, nsmul_eq_mul, nsmul_eq_mul] at hf
    push_cast at hf
    show f m / ((a : ℕ) : ℝ≥0) = f m' / ((a' : ℕ) : ℝ≥0)
    rw [div_eq_div_iff ha0 ha'0]
    refine mul_left_cancel₀ hk0 ?_
    calc ((k : ℕ) : ℝ≥0) * (f m * ((a' : ℕ) : ℝ≥0))
        = ((k : ℕ) : ℝ≥0) * ((a' : ℕ) : ℝ≥0) * f m := by ring
      _ = ((k : ℕ) : ℝ≥0) * ((a : ℕ) : ℝ≥0) * f m' := hf
      _ = ((k : ℕ) : ℝ≥0) * (f m' * ((a : ℕ) : ℝ≥0)) := by ring)

@[simp] theorem pfLiftFun_mk (f : M →+ ℝ≥0) (m : M) (a : ℕ+) :
    pfLiftFun f (Pf.mk m a) = f m / ((a : ℕ) : ℝ≥0) := rfl

/-- ★★**`⊗ ℝ≥0` の実体** —— `f : M →+ ℝ≥0` の `M^pf` への延長。

★`M^pf` の元は形式的な分数 `m/n` なので、`m/n ↦ f m / n` で延びる。 -/
noncomputable def pfLift (f : M →+ ℝ≥0) : Pf M →+ ℝ≥0 where
  toFun := pfLiftFun f
  map_zero' := by
    show pfLiftFun f (Pf.mk 0 1) = 0
    simp
  map_add' x y := by
    induction x using Pf.inductionOn with
    | _ m a =>
      induction y using Pf.inductionOn with
      | _ m' a' =>
        have ha0 : ((a : ℕ) : ℝ≥0) ≠ 0 := by
          simpa using (a.pos.ne' : (0 : ℕ) ≠ (a : ℕ))
        have ha'0 : ((a' : ℕ) : ℝ≥0) ≠ 0 := by
          simpa using (a'.pos.ne' : (0 : ℕ) ≠ (a' : ℕ))
        rw [Pf.mk_add_mk]
        show f ((a' : ℕ) • m + (a : ℕ) • m') / (((a * a' : ℕ+) : ℕ) : ℝ≥0)
          = f m / ((a : ℕ) : ℝ≥0) + f m' / ((a' : ℕ) : ℝ≥0)
        rw [map_add, map_nsmul, map_nsmul, nsmul_eq_mul, nsmul_eq_mul]
        push_cast
        field_simp

@[simp] theorem pfLift_mk (f : M →+ ℝ≥0) (m : M) (a : ℕ+) :
    pfLift f (Pf.mk m a) = f m / ((a : ℕ) : ℝ≥0) := rfl

@[simp] theorem pfLift_of (f : M →+ ℝ≥0) (m : M) : pfLift f (Pf.of m) = f m := by
  show pfLift f (Pf.mk m 1) = f m
  simp

/-- ★**`f` が単射なら延長も単射** —— `M^pf ↪ ℝ≥0`。

★`f m / a = f m' / a'` から `a' • m` と `a • m'` の `f` 値が一致し、
`f` の単射性で `a' • m = a • m'`、すなわち `Pf.mk m a = Pf.mk m' a'`。 -/
theorem pfLift_injective {f : M →+ ℝ≥0} (hf : Function.Injective f) :
    Function.Injective (pfLift f) := by
  intro x y hxy
  induction x using Pf.inductionOn with
  | _ m a =>
    induction y using Pf.inductionOn with
    | _ m' a' =>
      have ha0 : ((a : ℕ) : ℝ≥0) ≠ 0 := by
        simpa using (a.pos.ne' : (0 : ℕ) ≠ (a : ℕ))
      have ha'0 : ((a' : ℕ) : ℝ≥0) ≠ 0 := by
        simpa using (a'.pos.ne' : (0 : ℕ) ≠ (a' : ℕ))
      rw [pfLift_mk, pfLift_mk, div_eq_div_iff ha0 ha'0] at hxy
      refine Pf.sound 1 ?_
      refine hf ?_
      rw [map_nsmul, map_nsmul, nsmul_eq_mul, nsmul_eq_mul]
      push_cast
      calc (1 : ℝ≥0) * ((a' : ℕ) : ℝ≥0) * f m = f m * ((a' : ℕ) : ℝ≥0) := by ring
        _ = f m' * ((a : ℕ) : ℝ≥0) := hxy
        _ = (1 : ℝ≥0) * ((a : ℕ) : ℝ≥0) * f m' := by ring

/-! ## ★`Λ≥0 ↪ ℝ≥0` —— monoid type の実数化 -/

/-- ★`Λ≥0` から `ℝ≥0` への自然な埋め込み。

★`ℤ≥0 = ℕ ↪ ℝ≥0`、`ℚ≥0 ↪ ℝ≥0`、`ℝ≥0 = ℝ≥0`。 -/
noncomputable def MonoidType.toNNReal : (Λ : MonoidType) → (Λ.Nonneg →+ ℝ≥0)
  | .int => Nat.castAddMonoidHom ℝ≥0
  | .rat =>
      { toFun := fun q : ℚ≥0 => (q : ℝ≥0)
        map_zero' := NNRat.cast_zero
        map_add' := fun p q => NNRat.cast_add p q }
  | .real => AddMonoidHom.id ℝ≥0

theorem MonoidType.toNNReal_injective (Λ : MonoidType) :
    Function.Injective Λ.toNNReal := by
  cases Λ
  · exact fun a b h => Nat.cast_injective (R := ℝ≥0) h
  · exact fun a b h => NNRat.cast_injective h
  · exact fun a b h => h

/-- ★★**monoprime なモノイドの `M^pf` は `ℝ≥0` に埋まる**。

★★**これが `M^rlf_p = M^pf_p ⊗ ℝ≥0` の実体**である ——
`M ≃+ Λ≥0 ↪ ℝ≥0` を `pfLift` で `M^pf` へ延ばす。 -/
noncomputable def monoprimeLift {Λ : MonoidType} (e : M ≃+ Λ.Nonneg) : Pf M →+ ℝ≥0 :=
  pfLift ((Λ.toNNReal).comp (e : M →+ Λ.Nonneg))

theorem monoprimeLift_injective {Λ : MonoidType} (e : M ≃+ Λ.Nonneg) :
    Function.Injective (monoprimeLift e) :=
  pfLift_injective (fun _ _ h =>
    e.injective ((Λ.toNNReal_injective) h))

/-! ## ★★同型の取り方に依らないこと

★`monoprimeLift` は同型 `e` の選び方に依るが、★**依り方は `ℝ≥0` の正の実数倍だけ**
である。★これが原文が `⊗ ℝ≥0` を「標準的」と書ける根拠である。 -/

/-- ★★**延長は `M` 上の値で決まる** —— `M^pf` は `M` の像で生成されるので。 -/
theorem pfLift_ext {f g : M →+ ℝ≥0} (h : ∀ m : M, f m = g m) :
    pfLift f = pfLift g := by
  refine AddMonoidHom.ext fun x => ?_
  induction x using Pf.inductionOn with
  | _ m a => simp [h m]

/-- ★★★**正の実数倍を除く一意性** —— 2 つの延長が `M` 上で正の実数倍だけ違えば、
`M^pf` 上でも同じ倍率で違う。

★★**「`⊗ ℝ≥0` の取り方に依らない」の中身**である。 -/
theorem pfLift_smul {f g : M →+ ℝ≥0} (c : ℝ≥0) (h : ∀ m : M, g m = c * f m)
    (x : Pf M) : pfLift g x = c * pfLift f x := by
  induction x using Pf.inductionOn with
  | _ m a =>
    rw [pfLift_mk, pfLift_mk, h m, mul_div_assoc]

/-! ## ★段 1 の後半 —— `sup(Bound(a))`

原文 (FrdI p.47):
> the various Bound

★★**原文が `sup` と書ける根拠を式にする。**

★`sup` が取れるには **2 つ**要る:
1. **`M^rlf_p` が全順序かつ条件付き完備**であること
   —— ★**これが (b)(`M_p` は monoprime)の効き所**である。
   `M_p ≃ Λ≥0` から `(M_p)^pf` は `ℚ≥0` か `ℝ≥0`(`Λ = ℤ` なら `ℕ^pf ≃ ℚ≥0`)、
   したがって `⊗ ℝ≥0` して `M^rlf_p ≃ ℝ≥0` になる。
   ★**(b) が無ければ `M^pf_p ⊗ ℝ≥0` はただの加法モノイドで、`sup` は意味を持たない。**
2. その部分集合が **`M^rlf_p` の中で上に有界**であること
   —— ★★**これは自動ではない。** `Bound` は定義から `a` で押さえられているが、
   ★**その `a` は `M^pf` の元であって `M^pf_p` の元ではない**ので、
   `M^rlf_p` の中での上界にはならない。
   ★**だから原文は角括弧で「i.e., the various Bound(a) are bounded subsets」と
   条件として書いている。** ここでも `BddAbove` を仮定として持つ。
-/

/-- ★`ℝ≥0` ではモノイドの `≤`(`MLe`)は順序の `≤` に一致する。 -/
theorem mle_nnreal_iff {x y : ℝ≥0} : MLe x y ↔ x ≤ y := by
  constructor
  · rintro ⟨c, rfl⟩
    exact le_self_add
  · intro h
    exact ⟨y - x, add_tsub_cancel_of_le h⟩

/-- ★加法準同型は `MLe` を保つ。 -/
theorem map_mle {N : Type*} [AddCommMonoid N] (f : M →+ N) {a b : M} (h : MLe a b) :
    MLe (f a) (f b) := by
  obtain ⟨c, rfl⟩ := h
  exact ⟨f c, by rw [map_add]⟩

/-- ★`ℝ≥0` へ送れば `MLe` は `≤` になる。 -/
theorem map_le_of_mle (f : M →+ ℝ≥0) {a b : M} (h : MLe a b) : f a ≤ f b :=
  mle_nnreal_iff.mp (map_mle f h)

/-- ★★**(b) の効き所** —— monoprime なら `M^pf` は `ℝ≥0` へ
**順序を保って**単射に埋まる。★これで `M^rlf_p` は条件付き完備な全順序になる。 -/
theorem monoprimeLift_le_of_mle {Λ : MonoidType} (e : M ≃+ Λ.Nonneg) {a b : Pf M}
    (h : MLe a b) : monoprimeLift e a ≤ monoprimeLift e b :=
  map_le_of_mle _ h

variable (M) in
/-- ★`0` はつねに `Bound S a` に入る(`S` が `0` を含むなら)。

★★原文の添字が `p ∪ {0}` である理由 —— **`sup` を取る集合を空にしないため**。 -/
theorem zero_mem_bound {S : Set M} (h0 : (0 : M) ∈ S) (a : M) :
    (0 : M) ∈ Bound M S a :=
  ⟨h0, ⟨a, zero_add a⟩⟩

/-- ★★**原文の `sup(Bound_S(a))`** —— `Bound_S(a)` の像の上限を `ℝ≥0` で取る。 -/
noncomputable def boundSup (ι : M → ℝ≥0) (S : Set M) (a : M) : ℝ≥0 :=
  sSup (ι '' Bound M S a)

/-- ★上界であること。★**`BddAbove` が要る**(上の節の 2 番)。 -/
theorem le_boundSup (ι : M → ℝ≥0) (S : Set M) (a : M)
    (hbdd : BddAbove (ι '' Bound M S a)) {x : M} (hx : x ∈ Bound M S a) :
    ι x ≤ boundSup ι S a :=
  le_csSup hbdd ⟨x, hx, rfl⟩

/-- ★**最小の上界**であること。★`S` が `0` を含めば集合は空でない。 -/
theorem boundSup_le (ι : M → ℝ≥0) {S : Set M} (h0 : (0 : M) ∈ S) (a : M) {c : ℝ≥0}
    (h : ∀ x ∈ Bound M S a, ι x ≤ c) : boundSup ι S a ≤ c :=
  csSup_le (⟨ι 0, ⟨0, zero_mem_bound M h0 a, rfl⟩⟩) (by rintro _ ⟨x, hx, rfl⟩; exact h x hx)

/-- ★★`Bound` は単調 —— `MLe a b` なら `Bound S a ⊆ Bound S b`。 -/
theorem bound_mono {S : Set M} {a b : M} (h : MLe a b) : Bound M S a ⊆ Bound M S b :=
  fun _ hx => ⟨hx.1, mle_trans hx.2 h⟩

/-- ★★**`boundSup` は単調** —— 因子分解写像が `MLe` を保つことの中身。 -/
theorem boundSup_mono (ι : M → ℝ≥0) (S : Set M) {a b : M} (h : MLe a b)
    (hbdd : BddAbove (ι '' Bound M S b)) (h0 : (0 : M) ∈ S) :
    boundSup ι S a ≤ boundSup ι S b :=
  boundSup_le ι h0 a fun x hx => le_boundSup ι S b hbdd (bound_mono h hx)

/-! ## ★★段 2 —— `Definition 2.4, (i)` perf-factorial

原文 (FrdI p.47):
> Definition 2.4.

★★**族 `ι` をパラメータに出す**(案 B) —— 原文も `M^rlf_p` を
「(b) から得られているもの」として名前で参照し、条件の列挙とは別の層に置いている。
★段 1 の `pfLift_smul`(**正の実数倍を除いて一意**)がこの層を安全にしている。

★★**(c1)(c3) は「条件」であって定理ではない** ——
原文は (c) の中で
「well-defined [i.e., the various Bound(a) are bounded subsets]」
「whose image lies in ∏_p M^pf_p」
と**角括弧と関係詞で条件として**書いている。★我々の測定(2026-08-17)でも、
`Bound` を押さえる `a` は `M^pf` の元であって `M^pf_p` の元ではないので
有界性は自動ではなく、また `M^pf_p ≃ ℚ≥0` は `ℝ≥0` の中で `sup` について閉じない。
★**したがって条件として書くのが正しい。**
-/

/-- ★★`M^rlf_factor := ∏_{p ∈ Prime(M)} M^rlf_p`。

★段 1 により各 `M^rlf_p ≃ ℝ≥0` なので、`Prime(M) → ℝ≥0` として実現する。 -/
abbrev RlfFactor (M : Type w) [AddCommMonoid M] : Type w := Prime M → ℝ≥0

/-- ★**`Supp(a) ⊆ Prime(M)`** —— 成分が `0` でない素点の集合。 -/
def Supp (a : RlfFactor M) : Set (Prime M) := {p | a p ≠ 0}

variable (M) in
/-- ★**`M^pf` の中の「`p` ∪ {0}」** —— 原文 `Bound^p_{0}` の添字集合。

★`p` の元は「正の整数倍が `Primary(M)` の `p` 類に入る `M^pf` の元」
(§0 p.12 の `Primary(M^pf)` の特徴づけ、`isPrimaryElt_pf_iff`)。
★`{0}` を足すのは **`sup` を取る集合を空にしないため**。 -/
def pCarrierPf (p : Prime M) : Set (Pf M) :=
  {x : Pf M | ∃ (n : ℕ+) (b : M), b ∈ primeCarrier M p ∧ ((n : ℕ+) : ℕ) • x = Pf.of b}
    ∪ {0}

theorem zero_mem_pCarrierPf (p : Prime M) : (0 : Pf M) ∈ pCarrierPf M p :=
  Or.inr rfl

/-- ★★**因子分解写像**(族 `ι` を与えたとき)——
原文の `a ↦ (…, sup(Bound_{p ∪ {0}}(a)), …)`。 -/
noncomputable def factorMap (ι : Prime M → Pf M → ℝ≥0) (a : Pf M) : RlfFactor M :=
  fun p => boundSup (ι p) (pCarrierPf M p) a

/-- ★★★**[FrdI] Definition 2.4, (i)** —— `perf-factorial`(族 `ι` を明示した形)。

★**フィールドは原文の (a)(b)(c)(d) と 1 対 1**((c) は原文どおり 3 つに割る)。 -/
structure IsPerfFactorialWith (M : Type w) [AddCommMonoid M]
    (ι : Prime M → Pf M → ℝ≥0) : Prop where
  /-- **(a)** `M` is divisorial. -/
  divisorial : IsDivisorial M
  /-- **(b)** 各 `p ∈ Prime(M)` で `M_p` は monoprime。 -/
  monoprimeAt : ∀ p : Prime M, IsMonoprime (Mp M p)
  /-- ★`ι p` は `p` 成分の上で**加法的**(`M^rlf_p ≃ ℝ≥0` の同一視)。 -/
  embedAdd : ∀ (p : Prime M) {x y : Pf M}, x ∈ pCarrierPf M p → y ∈ pCarrierPf M p →
      ι p (x + y) = ι p x + ι p y
  /-- ★`ι p` は `p` 成分の上で**単射**。 -/
  embedInj : ∀ p : Prime M, Set.InjOn (ι p) (pCarrierPf M p)
  /-- ★`ι p` は `p` 成分の上で **`≤` を保つ**。 -/
  embedMono : ∀ (p : Prime M) {x y : Pf M}, x ∈ pCarrierPf M p → y ∈ pCarrierPf M p →
      MLe x y → ι p x ≤ ι p y
  /-- **(c1)** ★**原文が角括弧で条件として書く所** —— `Bound` の像が上に有界。 -/
  bounded : ∀ (a : Pf M) (p : Prime M),
      BddAbove (ι p '' Bound (Pf M) (pCarrierPf M p) a)
  /-- **(c2)** 因子分解写像は**加法的**。 -/
  factorAdd : ∀ a b : Pf M, factorMap ι (a + b) = factorMap ι a + factorMap ι b
  /-- **(c2)** 因子分解写像は**単射**。 -/
  factorInj : Function.Injective (factorMap ι)
  /-- **(c3)** ★**これも原文が条件として書く所** —— 像が `∏_p M^pf_p` に入る。 -/
  factorMem : ∀ (a : Pf M) (p : Prime M), factorMap ι a p ∈ ι p '' pCarrierPf M p
  /-- **(d)** `Supp` の条件 ——
  `a ∈ M^pf_factor`、`b ∈ M^pf`、`Supp(a) ⊆ Supp(b)` ならば `a ∈ M^pf`。 -/
  supp : ∀ a : RlfFactor M, (∀ p, a p ∈ ι p '' pCarrierPf M p) →
      ∀ b : Pf M, Supp a ⊆ Supp (factorMap ι b) → a ∈ Set.range (factorMap ι)

/-- ★★★**[FrdI] Definition 2.4, (i)** —— `perf-factorial`。

★族 `ι`(＝各素点での `M^rlf_p ≃ ℝ≥0` の同一視)を存在量化した形。 -/
def IsPerfFactorial (M : Type w) [AddCommMonoid M] : Prop :=
  ∃ ι : Prime M → Pf M → ℝ≥0, IsPerfFactorialWith M ι

/-- ★★**realification `M^rlf ⊆ M^rlf_factor`**(原文 p.48)——
`∃ b ∈ M^pf, Supp(a) ⊆ Supp(b)` を満たす `a` の全体。 -/
def Rlf (ι : Prime M → Pf M → ℝ≥0) : Set (RlfFactor M) :=
  {a : RlfFactor M | ∃ b : Pf M, Supp a ⊆ Supp (factorMap ι b)}

/-- ★`M^pf` の像は `M^rlf` に入る。 -/
theorem factorMap_mem_rlf (ι : Prime M → Pf M → ℝ≥0) (b : Pf M) :
    factorMap ι b ∈ Rlf ι :=
  ⟨b, subset_rfl⟩

/-! ### ★★族の取り替えに対する不変性

★★**案 B(族をパラメータに出す)を安全にしている根拠**である ——
段 1 の `pfLift_smul` により、族の取り方は**各素点で正の実数倍を除いて一意**。
★**その取り替えで `IsPerfFactorialWith` は保たれる。** -/

/-- ★`ℝ≥0` では正のスカラー倍は `sSup` と可換。 -/
theorem boundSup_smul (ι : M → ℝ≥0) (S : Set M) (a : M) {c : ℝ≥0} (hc : c ≠ 0)
    (h0 : (0 : M) ∈ S) (hbdd : BddAbove (ι '' Bound M S a)) :
    boundSup (fun x => c * ι x) S a = c * boundSup ι S a := by
  have hne : (ι '' Bound M S a).Nonempty := ⟨ι 0, 0, zero_mem_bound M h0 a, rfl⟩
  have h := (OrderIso.mulLeft₀ c (pos_iff_ne_zero.mpr hc)).map_csSup' hne hbdd
  show sSup ((fun x => c * ι x) '' Bound M S a) = c * sSup (ι '' Bound M S a)
  rw [show ((fun x => c * ι x) '' Bound M S a)
      = (OrderIso.mulLeft₀ c (pos_iff_ne_zero.mpr hc)) '' (ι '' Bound M S a) from by
    rw [Set.image_image]; rfl]
  exact h.symm

/-- ★因子分解写像は族のスカラー倍に**共変**。 -/
theorem factorMap_smul {ι ι' : Prime M → Pf M → ℝ≥0} (c : Prime M → ℝ≥0)
    (hc : ∀ p, c p ≠ 0) (h : ∀ p x, ι' p x = c p * ι p x)
    (hbdd : ∀ (a : Pf M) (p : Prime M),
      BddAbove (ι p '' Bound (Pf M) (pCarrierPf M p) a))
    (a : Pf M) (p : Prime M) :
    factorMap ι' a p = c p * factorMap ι a p := by
  show boundSup (ι' p) (pCarrierPf M p) a = c p * boundSup (ι p) (pCarrierPf M p) a
  rw [show ι' p = fun x => c p * ι p x from funext (h p)]
  exact boundSup_smul (ι p) _ a (hc p) (zero_mem_pCarrierPf p) (hbdd a p)

/-- ★★★**族の取り替えに対する不変性** ——
各素点で正の実数倍だけ違う族に取り替えても `perf-factorial` 性は保たれる。 -/
theorem IsPerfFactorialWith.smul {ι ι' : Prime M → Pf M → ℝ≥0} (c : Prime M → ℝ≥0)
    (hc : ∀ p, c p ≠ 0) (h : ∀ p x, ι' p x = c p * ι p x)
    (H : IsPerfFactorialWith M ι) : IsPerfFactorialWith M ι' := by
  have hfm : ∀ (a : Pf M) (p : Prime M), factorMap ι' a p = c p * factorMap ι a p :=
    factorMap_smul c hc h H.bounded
  have himg : ∀ (p : Prime M) (T : Set (Pf M)),
      ι' p '' T = (fun z => c p * z) '' (ι p '' T) := by
    intro p T
    rw [Set.image_image]
    exact congrArg (· '' T) (funext (h p))
  refine ⟨H.divisorial, H.monoprimeAt, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro p x y hx hy
    rw [h, h, h, H.embedAdd p hx hy, mul_add]
  · intro p x hx y hy hxy
    refine H.embedInj p hx hy ?_
    rw [h, h] at hxy
    exact mul_left_cancel₀ (hc p) hxy
  · intro p x y hx hy hle
    rw [h, h]
    exact mul_le_mul_left' (H.embedMono p hx hy hle) _
  · intro a p
    obtain ⟨u, hu⟩ := H.bounded a p
    refine ⟨c p * u, ?_⟩
    rintro _ ⟨x, hx, rfl⟩
    rw [h]
    exact mul_le_mul_left' (hu ⟨x, hx, rfl⟩) _
  · intro a b
    refine funext fun p => ?_
    simp only [Pi.add_apply]
    rw [hfm, hfm, hfm, congrFun (H.factorAdd a b) p]
    simp only [Pi.add_apply]
    exact mul_add _ _ _
  · intro a b hab
    refine H.factorInj (funext fun p => ?_)
    have hp := congrFun hab p
    rw [hfm, hfm] at hp
    exact mul_left_cancel₀ (hc p) hp
  · intro a p
    rw [hfm, himg p (pCarrierPf M p)]
    exact ⟨factorMap ι a p, H.factorMem a p, rfl⟩
  · intro a ha b hsupp
    have ha' : ∀ p, ((c p)⁻¹ * a p) ∈ ι p '' pCarrierPf M p := by
      intro p
      obtain ⟨x, hx, hxe⟩ := ha p
      refine ⟨x, hx, ?_⟩
      rw [h p x] at hxe
      rw [← hxe, ← mul_assoc, inv_mul_cancel₀ (hc p), one_mul]
    have hsupp' : Supp (fun q => (c q)⁻¹ * a q) ⊆ Supp (factorMap ι b) := by
      intro p hp
      have hap : a p ≠ 0 := by
        intro h0
        exact hp (show (c p)⁻¹ * a p = 0 by rw [h0, mul_zero])
      have hin : p ∈ Supp (factorMap ι' b) := hsupp hap
      show factorMap ι b p ≠ 0
      intro h0
      exact hin (by rw [hfm, h0, mul_zero])
    obtain ⟨y, hy⟩ := H.supp _ ha' b hsupp'
    refine ⟨y, funext fun p => ?_⟩
    rw [hfm, congrFun hy p, ← mul_assoc, mul_inv_cancel₀ (hc p), one_mul]

/-! ## ★段 3 —— `Definition 2.4, (ii)` `Λ supports M`

原文 (FrdI p.48):
> (ii) Let

★**3 つの場合分けそのもの**である。 -/

/-- **[FrdI] Definition 2.4, (ii)** —— `Λ>0`。 -/
def MonoidType.Pos : MonoidType → Type
  | .int => ℕ+
  | .rat => {q : ℚ≥0 // q ≠ 0}
  | .real => {r : ℝ≥0 // r ≠ 0}

/-- ★★★**[FrdI] Definition 2.4, (ii)** —— **`Λ supports M`**。

★`(a) Λ = ℤ`、`(b) Λ = ℚ` かつ `M` perfect、
`(c) Λ = ℝ` かつ `M` perfect かつ perf-factorial かつ各 `M_p` が ℝ-monoprime。 -/
def Supports (Λ : MonoidType) (M : Type w) [AddCommMonoid M] : Prop :=
  match Λ with
  | .int => True
  | .rat => IsPerfectMonoid M
  | .real => IsPerfectMonoid M ∧ IsPerfFactorial M
      ∧ ∀ p : Prime M, IsLambdaMonoprime (Mp M p) MonoidType.real

/-- ★`ℤ` はつねに supports する。 -/
theorem supports_int (M : Type w) [AddCommMonoid M] : Supports MonoidType.int M := trivial

/-- ★`n •` を加法準同型として。 -/
def nsmulHom (n : ℕ) (M : Type w) [AddCommMonoid M] : M →+ M where
  toFun x := n • x
  map_zero' := smul_zero n
  map_add' x y := smul_add n x y

@[simp] theorem nsmulHom_apply (n : ℕ) (M : Type w) [AddCommMonoid M] (x : M) :
    nsmulHom n M x = n • x := rfl

/-- ★★**`M` が perfect なら `n •` は加法同型** —— これが `Λ = ℚ` の作用の一意性の中身。 -/
noncomputable def nsmulEquiv {M : Type w} [AddCommMonoid M] (h : IsPerfectMonoid M)
    (n : ℕ+) : M ≃+ M :=
  AddEquiv.ofBijective (nsmulHom ((n : ℕ+) : ℕ) M) (h n)

/-! ### ★「`Λ>0` が `M` に自然に作用する」について

★原文は (ii) の末尾で
「Note that if `Λ` supports `M`, then `Λ>0` acts naturally on `M`」
と注意する。★**その作用の中身**:

- **`Λ = ℤ`**: `n • x` そのもの(`nsmulHom`)。
- **`Λ = ℚ`**: `q = c/e` に対し `q · x` は `e • y = c • x` なる一意な `y`。
  ★**一意性を与えるのが `M` の perfect 性**(`nsmulEquiv`)。
- **`Λ = ℝ`**: 因子分解 `M^pf ↪ ∏_p M^rlf_p` の各成分で `ℝ≥0` のスカラー倍を取る。
  ★★**ここで「各 `M_p` が ℝ-monoprime」が効く** —— 成分が `ℝ≥0` **全体**なので
  正のスカラー倍で閉じ、`Supp` が変わらないので **(d) の条件により `M^pf` に戻る**。
  ★**(c) が ℝ-monoprime を要求する理由がこれである。**
-/

/-! ## ★段 4 —— `Definition 2.4, (iii)` `d·Φ(−)` / `𝒞(d)` / `d` の Frobenius 函手

原文 (FrdI p.48):
> for the subcategory determined by the arrows whose zero divisor lies in d

★**以下は `Λ = ℤ`(すなわち `d ∈ ℕ≥1`)の場合を構成する** ——
`Proposition 2.1, (ii)` と同じ設定であり、`Proposition 2.5, (iii)` が使うのもこの形。 -/

section Cd


variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)


variable {P} in
/-- ★★**`d · Φ(−) ⊆ Φ(−)`**(`Definition 2.4, (iii)`)。 -/
def scalePhi (d : ℕ+) (Y : D) : AddSubmonoid (Φ.val Y) :=
  AddMonoidHom.mrange (nsmulHom ((d : ℕ+) : ℕ) (Φ.val Y))

theorem mem_scalePhi_iff (d : ℕ+) {Y : D} (x : Φ.val Y) :
    x ∈ scalePhi (Φ := Φ) d Y ↔ ∃ u, ((d : ℕ+) : ℕ) • u = x := Iff.rfl

/-- ★**部分関手であること** —— `Φ.map` は `d · Φ(−)` を保つ。 -/
theorem scalePhi_map (d : ℕ+) {Y Z : D} (α : Z ⟶ Y) {x : Φ.val Y}
    (hx : x ∈ scalePhi (Φ := Φ) d Y) : Φ.map α x ∈ scalePhi (Φ := Φ) d Z := by
  obtain ⟨u, rfl⟩ := hx
  exact ⟨Φ.map α u, by rw [nsmulHom_apply, nsmulHom_apply, map_nsmul]⟩

/-- ★★**`𝒞(d)` を定める射の性質** —— 零因子が `d · Φ(−)` に入ること。 -/
def cdProp (d : ℕ+) : MorphismProperty C :=
  fun A _ φ => P.Div φ ∈ scalePhi (Φ := Φ) d (P.toElem.obj A).base

instance cdProp_isMultiplicative (d : ℕ+) : (cdProp P d).IsMultiplicative where
  id_mem A := ⟨0, by rw [nsmulHom_apply, smul_zero]; exact (P.Div_id A).symm⟩
  comp_mem {A B E} ψ φ hψ hφ := by
    obtain ⟨v, hv⟩ := hψ
    obtain ⟨u, hu⟩ := hφ
    refine ⟨Φ.map (P.Base ψ) u + ((P.degFr φ : ℕ+) : ℕ) • v, ?_⟩
    rw [nsmulHom_apply] at hv hu ⊢
    rw [P.Div_comp, ← hv, ← hu, smul_add, map_nsmul, smul_comm]

/-- ★★★**`𝒞(d) ⊆ 𝒞`**(`Definition 2.4, (iii)`)。 -/
abbrev Cd (d : ℕ+) : Type u2 := WideSubcategory (cdProp P d)

/-- ★★**`Φ` の自己準同型「`d` 倍」** —— `Definition 1.1, (iii)` の意味での
`𝔽_Φ → 𝔽_Φ` を誘導する自然変換。 -/
def phiNsmulNat (d : ℕ+) : Φ.functor ⟶ Φ.functor where
  app Y := AddCommMonCat.ofHom (nsmulHom ((d : ℕ+) : ℕ) _)
  naturality {Y Z} f := by
    refine AddCommMonCat.hom_ext (AddMonoidHom.ext fun x => ?_)
    show ((d : ℕ+) : ℕ) • (Φ.functor.map f).hom x
      = (Φ.functor.map f).hom (((d : ℕ+) : ℕ) • x)
    rw [map_nsmul]

/-- ★★★**`d` の Frobenius 函手 `𝔽_Φ → 𝔽_Φ`**(`Definition 2.4, (iii)`)。

★**`Proposition 2.1, (ii)` の `𝔽_Φ` 側の函手そのもの**である。 -/
def frobFunctorOfDeg (d : ℕ+) : ElemFrobCat Φ ⥤ ElemFrobCat Φ :=
  ElemFrobCat.elemFrobMap (phiNsmulNat (Φ := Φ) d)

@[simp] theorem frobFunctorOfDeg_obj (d : ℕ+) (A : ElemFrobCat Φ) :
    (frobFunctorOfDeg (Φ := Φ) d).obj A = ⟨A.base⟩ := rfl

/-- ★**Frobenius 次数と両立する**(原文の "compatible with Frobenius degrees")。 -/
@[simp] theorem frobFunctorOfDeg_deg (d : ℕ+) {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    ((frobFunctorOfDeg (Φ := Φ) d).map φ).deg = φ.deg := rfl

/-- ★**`𝒟` への射影と両立する**(原文の "the natural projection functor")。 -/
@[simp] theorem frobFunctorOfDeg_base (d : ℕ+) {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    ((frobFunctorOfDeg (Φ := Φ) d).map φ).base = φ.base := rfl

@[simp] theorem frobFunctorOfDeg_div (d : ℕ+) {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    ((frobFunctorOfDeg (Φ := Φ) d).map φ).div = ((d : ℕ+) : ℕ) • φ.div := rfl

/-- ★★**像が `(𝔽_Φ)(d)` に入る** —— 原文の `C(d) → F_{d·Φ} = (F_Φ)(d)`。 -/
theorem frobFunctorOfDeg_div_mem (d : ℕ+) {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    ((frobFunctorOfDeg (Φ := Φ) d).map φ).div ∈ scalePhi (Φ := Φ) d A.base :=
  ⟨φ.div, rfl⟩

/-! ### ★★一般の `Λ` —— 作用をパラメータにした形

★原文の (iii) は `d ∈ Λ>0` について述べる。★**`Λ>0` の作用**(＝ (ii) の末尾の注意)
**を与えれば、`d · Φ(−)` と `𝒞(d)` は同じ形で書ける。**
★`Λ = ℤ` の場合が上の `scalePhi` / `cdProp` である。 -/

variable {P} in
/-- ★★`Φ` 上のスカラー作用(`d ∈ Λ>0` が定める自己準同型の族)。

★**要求は自然性だけ**である —— `ℕ` 倍との可換性は加法準同型から自動で従う。 -/
structure IsScalarAction (σ : ∀ Y : D, Φ.val Y →+ Φ.val Y) : Prop where
  /-- ★`Φ.map` と可換すること(部分関手になるための条件)。 -/
  natural : ∀ {Y Z : D} (α : Z ⟶ Y) (x : Φ.val Y), Φ.map α (σ Y x) = σ Z (Φ.map α x)

variable {P} in
/-- ★★**`d · Φ(−) ⊆ Φ(−)`**(一般の `Λ`)。 -/
def scalePhiOf (σ : ∀ Y : D, Φ.val Y →+ Φ.val Y) (Y : D) : AddSubmonoid (Φ.val Y) :=
  AddMonoidHom.mrange (σ Y)

/-- ★**部分関手であること**。 -/
theorem scalePhiOf_map {σ : ∀ Y : D, Φ.val Y →+ Φ.val Y} (hσ : IsScalarAction σ)
    {Y Z : D} (α : Z ⟶ Y) {x : Φ.val Y} (hx : x ∈ scalePhiOf σ Y) :
    Φ.map α x ∈ scalePhiOf σ Z := by
  obtain ⟨u, rfl⟩ := hx
  exact ⟨Φ.map α u, (hσ.natural α u).symm⟩

/-- ★★**`𝒞(d)` を定める射の性質**(一般の `Λ`)。 -/
def cdPropOf (σ : ∀ Y : D, Φ.val Y →+ Φ.val Y) : MorphismProperty C :=
  fun A _ φ => P.Div φ ∈ scalePhiOf σ (P.toElem.obj A).base

theorem cdPropOf_isMultiplicative {σ : ∀ Y : D, Φ.val Y →+ Φ.val Y}
    (hσ : IsScalarAction σ) : (cdPropOf P σ).IsMultiplicative where
  id_mem A := ⟨0, by rw [map_zero]; exact (P.Div_id A).symm⟩
  comp_mem {A B E} ψ φ hψ hφ := by
    obtain ⟨v, hv⟩ := hψ
    obtain ⟨u, hu⟩ := hφ
    refine ⟨Φ.map (P.Base ψ) u + ((P.degFr φ : ℕ+) : ℕ) • v, ?_⟩
    rw [map_add, map_nsmul, ← hσ.natural, hu, hv, P.Div_comp]

/-- ★★★**`𝒞(d) ⊆ 𝒞`**(一般の `Λ`)。 -/
abbrev CdOf (σ : ∀ Y : D, Φ.val Y →+ Φ.val Y) (hσ : IsScalarAction σ) : Type u2 :=
  letI := cdPropOf_isMultiplicative P hσ
  WideSubcategory (cdPropOf P σ)

/-- ★`Λ = ℤ` の作用(`d •`)は確かにスカラー作用である。 -/
theorem isScalarAction_nsmul (d : ℕ+) :
    IsScalarAction (Φ := Φ) (fun Y => nsmulHom ((d : ℕ+) : ℕ) (Φ.val Y)) :=
  ⟨fun α x => by rw [nsmulHom_apply, nsmulHom_apply, map_nsmul]⟩

end Cd

end ABC3.Found.FrdI
