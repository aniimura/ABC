/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24RlfPerf
import ABC3.Found.FrdI.Thm52

/-!
# `(M^pf)^gp ≅ (M^gp)^pf` —— 完全化と群化は交換する

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.104。

原文 (FrdI p.104):
> then so is Cpf. Moreover, Cun-tr, Crlf are always of model type. Finally, suppose

## ★なぜ要るか

`Proposition 5.5, (iv)` は「`𝒞` が model Frobenioid `(Φ, B)` なら
`𝒞^pf` は model Frobenioid `(Φ^pf, B^pf)`」と言う。
★model Frobenioid の定義は**有理関数の単系 `B` が `Φ^gp` の部分群**であることを
要求するので、`B^pf` を `(Φ^pf)^gp` の部分群として読むには

  `(Φ^gp)^pf ≅ (Φ^pf)^gp`

が要る(依存グラフの節点 `p55iv-pf` の「`B^pf → (Φ^gp)^pf` と `(Φ^pf)^gp` の比較」)。

## ★★逆写像を作らずに済ませる

★全単射を**直接**示す:

* **全射** —— `M^gp` の元は差 `toGp m − toGp m'` なので(`gp_eq_sub'`)、
  `mk z a` の原像は `toGp (m/a) − toGp (m'/a)` である。
* **単射** —— `(M^pf)^gp` の元も差 `toGp p − toGp q` であり、
  像が `0` なら `Pf.map` の単射性(`Pf.map_injective`、`M` が integral)で `p = q`。

★★これで**逆写像の well-definedness を書かずに済む**
(`Gp (Pf M)` の可除性を経由する必要がない)。

## ★本ファイルで閉じること

| 宣言 | 中身 |
|---|---|
| `Pf.negFun` / `Pf.instAddCommGroup` | ★**群の完全化は群** |
| `Pf.mk_sub_mk` | `g/a − g'/a = (g−g')/a` |
| `gp_eq_sub'` | `M^gp` の元は差 |
| `pfGpHom` | ★**`(M^pf)^gp ⟶ (M^gp)^pf`** |
| `pfGpEquiv` | ★★★**同型**(`M` が integral のとき) |
-/

namespace ABC3.Found.FrdI

universe w

/-! ## ★1. 群の完全化は群 -/

namespace Pf

variable {G : Type w} [AddCommGroup G]

/-- ★分子を否定するだけ。 -/
def negFun : Pf G → Pf G :=
  Quotient.map (fun x : G × ℕ+ => (-x.1, x.2)) (by
    rintro ⟨g, a⟩ ⟨g', a'⟩ ⟨k, e⟩
    refine ⟨k, ?_⟩
    have := congrArg (fun t : G => -t) e
    simpa using this)

instance : Neg (Pf G) := ⟨negFun⟩

@[simp] theorem neg_mk (g : G) (a : ℕ+) : -(mk g a) = mk (-g) a := rfl

noncomputable instance instAddCommGroup : AddCommGroup (Pf G) where
  __ := (inferInstance : AddCommMonoid (Pf G))
  neg := negFun
  zsmul := zsmulRec
  neg_add_cancel x := by
    induction x using Pf.inductionOn with | _ g a =>
    show mk (-g) a + mk g a = 0
    rw [mk_add_mk]
    refine Pf.sound 1 ?_
    simp

/-- ★同じ分母どうしの差。 -/
theorem mk_sub_mk (g g' : G) (a : ℕ+) : mk g a - mk g' a = mk (g - g') a := by
  rw [sub_eq_add_neg, neg_mk, mk_add_mk]
  refine Pf.sound 1 ?_
  push_cast
  simp only [smul_neg, smul_add, smul_smul, one_mul, smul_sub]
  abel

end Pf

/-! ## ★2. `M^gp` の元は差 -/

variable {M : Type w} [AddCommMonoid M]

/-- ★**`M^gp` の元はつねに `toGp a − toGp b` の形**。

★`mk_add_toGp`(群での消去)そのものである。 -/
theorem gp_eq_sub' (z : Gp M) : ∃ a b : M, z = toGp M a - toGp M b := by
  induction z using AddLocalization.induction_on with | _ x =>
  exact ⟨x.1, (x.2 : M), by rw [eq_sub_iff_add_eq]; exact mk_add_toGp M x.1 x.2⟩

/-- ★`M` が integral ⟺ `toGpHom` が単射。 -/
theorem toGpHom_injective (h : IsIntegralMonoid M) : Function.Injective (toGpHom M) := h

/-! ## ★3. `(M^pf)^gp ⟶ (M^gp)^pf` -/

/-- ★★**比較写像 `(M^pf)^gp ⟶ (M^gp)^pf`** ——
`M^pf ⟶ (M^gp)^pf`(`Pf.map (toGpHom M)`)を群化の普遍性で伸ばしたもの。 -/
noncomputable def pfGpHom : Gp (Pf M) →+ Pf (Gp M) :=
  gpLiftTo (Pf.map (toGpHom M))

@[simp] theorem pfGpHom_toGp (p : Pf M) :
    pfGpHom (toGp (Pf M) p) = Pf.map (toGpHom M) p :=
  gpLiftTo_toGp _ p

/-- ★★**全射** —— `mk z a` の原像は `toGp (m/a) − toGp (m'/a)`。 -/
theorem pfGpHom_surjective : Function.Surjective (pfGpHom (M := M)) := by
  intro y
  induction y using Pf.inductionOn with | _ z a =>
  obtain ⟨m, m', rfl⟩ := gp_eq_sub' z
  refine ⟨toGp (Pf M) (Pf.mk m a) - toGp (Pf M) (Pf.mk m' a), ?_⟩
  rw [map_sub, pfGpHom_toGp, pfGpHom_toGp, Pf.map_mk, Pf.map_mk, Pf.mk_sub_mk]
  rfl

/-- ★★**単射**(`M` が integral のとき)—— `Pf.map` の単射性 1 本。 -/
theorem pfGpHom_injective (h : IsIntegralMonoid M) :
    Function.Injective (pfGpHom (M := M)) := by
  refine (injective_iff_map_eq_zero _).mpr fun y hy => ?_
  obtain ⟨p, q, rfl⟩ := gp_eq_sub' y
  rw [map_sub, pfGpHom_toGp, pfGpHom_toGp, sub_eq_zero] at hy
  rw [Pf.map_injective (toGpHom_injective h) hy, sub_self]

/-- ★★★★★**完全化と群化は交換する** —— `(M^pf)^gp ≅ (M^gp)^pf`。

★これで `Proposition 5.5, (iv)` の `B^pf ⊆ (Φ^pf)^gp` が読める。 -/
noncomputable def pfGpEquiv (h : IsIntegralMonoid M) : Gp (Pf M) ≃+ Pf (Gp M) :=
  AddEquiv.ofBijective pfGpHom ⟨pfGpHom_injective h, pfGpHom_surjective⟩

@[simp] theorem pfGpEquiv_toGp (h : IsIntegralMonoid M) (p : Pf M) :
    pfGpEquiv h (toGp (Pf M) p) = Pf.map (toGpHom M) p :=
  pfGpHom_toGp p

/-! ## ★4. 自然性 -/

/-- ★`M^gp` からの射は `toGp` 上の値で決まる。 -/
theorem gp_hom_ext {G : Type*} [AddCommGroup G] {f g : Gp M →+ G}
    (h : ∀ m : M, f (toGp M m) = g (toGp M m)) : f = g := by
  ext z
  obtain ⟨a, b, rfl⟩ := gp_eq_sub' z
  rw [map_sub, map_sub, h a, h b]

/-- ★★★**比較写像は自然** —— `f : M →+ N` に沿って四角が可換。

★★これが `Φ` 上の同型(`MonoidOn` の同型)へ持ち上げるときに要る。 -/
theorem pfGpHom_naturality {N : Type w} [AddCommMonoid N] (f : M →+ N) :
    (pfGpHom (M := N)).comp (gpMap (Pf M) (Pf.map f))
      = (Pf.map (gpMap M f)).comp (pfGpHom (M := M)) := by
  refine gp_hom_ext fun p => ?_
  induction p using Pf.inductionOn with
  | _ m a =>
    show pfGpHom (gpMap (Pf M) (Pf.map f) (toGp (Pf M) (Pf.mk m a)))
      = Pf.map (gpMap M f) (pfGpHom (toGp (Pf M) (Pf.mk m a)))
    rw [gpMap_toGp, pfGpHom_toGp, pfGpHom_toGp]
    simp only [Pf.map_mk, toGpHom_apply, gpMap_toGp]

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.5, (iv)` が要する `(Φ^gp)^pf ≅ (Φ^pf)^gp`。 -/
def pfGpEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (iv) — (Φ^gp)^pf ≅ (Φ^pf)^gp(B^pf を読むため)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
